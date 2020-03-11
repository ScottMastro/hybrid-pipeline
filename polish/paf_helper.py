import pandas as pd
import variants as var
import external_tools as tools
import helper

def parse_paf(pafFile):
    """Loads file containing minimap2 alignments."""

    colnames = ["qid", "qlen", "qstart", "qend",
                "strand", 
                "rid", "rlen", "rstart", "rend", 
                "nmatch", "alnlen", "qual"]
    
    colorder = [c for c in colnames]
    
    typeMap = {"qlen": int, "qstart": int, "qend": int,
                    "rlen": int, "rstart": int, "rend": int,
                    "nmatch": int, "alnlen": int, "qual": int}
    
    lines = []
    with open(pafFile) as paf:
       line = paf.readline()
       while line:
           cols = line.split("\t")
           l = {c : x for c,x in zip(colnames, cols)}
           
           for i in range(len(colnames), len(cols)):
               col = cols[i]
               l[col[0:2]] = col[5:].rstrip("\n\r")
               if col[0:2] not in colorder: colorder.append(col[0:2])
               
           lines.append(l)
           line = paf.readline()
    
    
    df = pd.DataFrame(lines)
    df = df[colorder]
    
    if "dv" in df.columns:
        typeMap["dv"] = float    
    if "NM" in df.columns:
        typeMap["NM"] = int    


    df = df.astype(typeMap)

    return df


def parse_cs_string(cs, chrom, start, seq, noComplex=False): 
    '''The cs SAM/PAF tag encodes bases at mismatches and INDELs.
    It matches regular expression /(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/.
    Like CIGAR, cs consists of series of operations. 
    Each leading character specifies the operation; 
    the following sequence is the one involved in the operation.
    
    CGATCGATAAATAGAGTAG---GAATAGCA
    ||||||   ||||||||||   |||| |||
    CGATCG---AATAGAGTAGGTCGAATtGCA
    
    is represented as :6-ata:10+gtc:4*at:3, where :[0-9]+ represents an identical block,
    -ata represents a deletion, +gtc an insertion and *at indicates reference base a
    is substituted with a query base t. 
    It is similar to the MD SAM tag but is standalone and easier to parse.
    
    If --cs=long is used, the cs string also contains identical sequences in the alignment.
    The above example will become =CGATCG-ata=AATAGAGTAG+gtc=GAAT*at=GCA.
    The long form of cs encodes both reference and query sequences in one string.
    The cs tag also encodes intron positions and splicing signals (see the minimap2 manpage for details).
    '''
    
    variantSet = var.VariantCallSet()
    
    def make_variant_call(chrom, pos, ref, alt, annotation):
        allele = var.Allele(chrom, pos, ref, alt)
        return var.VariantCall(allele)

    def add_call(newVariantCall):
        #make call compound if a call already exists at the position we are adding to
        pos = newVariantCall.pos
        existing = variantSet.pop(newVariantCall.pos)
        if existing is not None:
            existingAllele = existing.main_allele()
            newAllele = newVariantCall.main_allele()
            chrom = existingAllele.chrom
            ref = existingAllele.ref + newAllele.ref[1:]
            alt = existingAllele.alt + newAllele.alt[1:]
            newVariantCall = make_variant_call(chrom, pos, ref, alt, "")
        variantSet.add(newVariantCall)
        
    EQL, INS, DEL, SUB = ":", "+", "-", "*"
    ops = [EQL, INS, DEL, SUB]
    
    currentPos = start
    
    startIdx, op = None, None
    pos, ref, alt, opList = None, "", "", ""


    for i,char in enumerate(cs):
        if char in ops:
            if op is not None:
                
                #to avoid outputting complex variants
                if noComplex:
                    # ...except if there is a gap being represented
                    if 'N' not in ref and 'N' not in alt:
                        if pos is not None:
                            add_call(make_variant_call(chrom, pos, ref, alt, opList))
                        pos, ref, alt, opList = None, "", "", ""
                
                val = cs[startIdx:i].upper()
                if op == EQL:
                    currentPos += int(val)
                    if pos is not None:
                        add_call(make_variant_call(chrom, pos, ref, alt, opList))
                    pos, ref, alt, opList = None, "", "", ""
                elif op == INS:
                    opList = opList + op
                    if pos is None:
                        anchor = seq[currentPos - 1]
                        pos = currentPos - 1
                        ref = anchor
                        alt = anchor + val
                    else:
                        alt = alt + val
                    #currentPos += 0
                elif op == DEL:
                    opList = opList + op
                    if pos is None:        
                        anchor = seq[currentPos - 1]
                        pos = currentPos - 1
                        ref = anchor + val
                        alt = anchor
                    else:
                        ref = ref + val
                    currentPos += len(val)
                elif op == SUB:
                    opList = opList + op
                    if pos is None:
                        pos = currentPos
                        ref = val[0]
                        alt = val[1]
                    else:
                        ref = ref + val[0]
                        alt = alt + val[1]
                    currentPos += 1
                        
            op = char
            startIdx=i+1
        
    if pos is not None and op != EQL:
        add_call(make_variant_call(chrom, pos, ref, alt, opList))
    
    return (variantSet)

def get_alignment_bounds(refFa, otherFa, prefix):
        
    pafFile = tools.align_paf_lenient(refFa, otherFa, prefix)
    pafDf = parse_paf(pafFile)
    #expectedSize = helper.get_fasta_len(otherFa)
    
    if len(pafDf) > 1:
        print("WARNING: multiple SV alignments!!")
        input()

        # could be a large insertion in the middle that is falsely aligned elsewhere
        # we trust the ends to be accurate
        firstAln = pafDf.sort_values(by=["qstart"], ascending=True).iloc[0]
        lastAln = pafDf.sort_values(by=["qend"], ascending=False).iloc[0]
        starts = [firstAln["qstart"], lastAln["qstart"]]
        otherAln = pafDf[~pafDf["qstart"].isin(starts)]
        
        # de = Gap-compressed per-base sequence divergence
        endsDv = [firstAln["de"], lastAln["de"]]
        middleDv = [aln["de"] for _,aln in otherAln.iterrows()]
        
        print("ends:", endsDv)
        print("middle:", middleDv)

        # Exclude middle alignments:
        pafDf = pafDf[pafDf["qstart"].isin(starts)]      

    qstart = min(min(pafDf["qstart"]), min(pafDf["qend"]))
    qend = max(max(pafDf["qstart"]), max(pafDf["qend"]))
    rstart = min(min(pafDf["rstart"]), min(pafDf["rend"]))
    rend = max(max(pafDf["rstart"]), max(pafDf["rend"]))

    return ((rstart,rend), (qstart,qend))


def get_variantset(pafFile, refFa, altFa):
    
    pafDf = parse_paf(pafFile)
    pafDf = pafDf.sort_values(by=["qstart"], ascending=True)
    
    refFaDict = tools.fasta2dict(refFa)    
    altFaDict = tools.fasta2dict(altFa)
        
    variantSetAll = var.VariantCallSet()
    qpos, rpos = None, None

    for _,aln in pafDf.iterrows():
        
        rid, qid = aln["rid"], aln["qid"]
        refSeq = refFaDict[rid]
        
        if qpos is not None:
            #todo: strand direction might mess this up
            #todo: separate chromosome alignment WILL mess this up

            qposEnd = aln["qstart"]
            rposEnd = aln["rstart"]

            ref = altFaDict[qid][qpos:qposEnd]
            alt = refSeq[rpos:rposEnd]

            gapCall = var.create_variant_call(rid, rpos, ref, alt)           
            variantSetAll.add(gapCall)
            
        aln["cs"]
        
        variantSet = parse_cs_string(aln["cs"], rid, int(aln["rstart"]), refSeq, noComplex=True)
        variantSetAll.add_all(variantSet)

        qpos, rpos = aln["qend"], aln["rend"]
        
    return variantSetAll

def get_variantset2(refFa, otherFa, prefix):
    
    pafFile = tools.align_paf_lenient(refFa, otherFa, prefix)
    pafDf = parse_paf(pafFile)
    pafDf = pafDf.sort_values(by=["qstart"], ascending=True)
    
    rid = pafDf.iloc[0]["rid"]
    faDict = tools.fasta2dict(refFa)
    seq = faDict[rid]

    firstAln = pafDf.sort_values(by=["qstart"], ascending=True).iloc[0]
    lastAln = pafDf.sort_values(by=["qend"], ascending=False).iloc[0]
    vs = parse_cs_string(firstAln["cs"], mainContig, int(firstAln["rstart"]), seqData[mainContig], noComplex=True)
    
    if len(pafDf) == 1:
        return vs
    
    vs2 = paf.parse_cs_string(lastAln["cs"], mainContig, int(lastAln["rstart"]), seqData[mainContig], noComplex=True)
    vs.add_all(vs2)
    
    if firstAln["qend"] > lastAln["qstart"]:
        print("WARNING: overlapping alignment - 1!!")
        #todo:fix
        input()
    
    if firstAln["rend"] > lastAln["rstart"]:
        print("WARNING: overlapping alignment - 2!!")
        #todo:fix
        input()

    altSeq = helper.get_fasta_seq(fa)
    seqRef = seqData[mainContig][firstAln["rend"]-1:lastAln["rstart"]]
    seqAlt = altSeq[firstAln["qend"]-1:lastAln["qstart"]]
    
    call = var.create_variant_call(mainContig, firstAln["rend"]-1, seqRef, seqAlt)
    vs.add(call)
    return vs


    #mainPaf = tools.align_paf_lenient(tigFa, mainPolishedFa, param.OUTPUT_DIR + "/" + mainContig + ".main")
    svPaf = tools.align_paf_lenient(mainPolishedFa, svPolishedFa, param.OUTPUT_DIR + "/" + mainContig + ".sv")
    
    #mainPafDf = paf.parse_paf(mainPaf)
    svPafDf = paf.parse_paf(svPaf)
    

    
    #variantSetMain = get_variantset(mainPafDf, mainPolishedFa)
    variantSetSV = get_variantset(svPafDf, svPolishedFa)
    #variantSet = var.combine_as_genotypes(variantSetMain, variantSetSV)
        
    '''
    outputAssessment = True
    if outputAssessment:
        fasta = helper.copy_file(regionFa, prefix + "_ASSESS_ref.fa")
        tools.align_pacbio(fasta, mainBam, prefix + "_ASSESS_mainreads")
        tools.align_pacbio(fasta, svBam, prefix + "_ASSESS_svreads")
        tools.align_pacbio(fasta, bamUnphased, prefix + "_ASSESS_unphased")

        tools.align_pacbio(fasta, mainPolishedFa,  prefix + "_ASSESS_maintig")
        tools.align_pacbio(fasta, svPolishedFa, prefix + "_ASSESS_svtig")
        variantSet = var.combine_as_genotypes(variantSetMain, variantSetSV)
        variantSet.write_vcf(prefix + "_ASSESS_variantset")
    ''' 
    
    #result.add_variants(var.VariantCallSet(), variantSetSV)
