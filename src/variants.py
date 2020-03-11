import subprocess
from enum import Enum
import copy

class Result(Enum):
 
    AGREE=1
    #agree but one or more is unphased
    AGREE_UNPHASED=2
    #agree on genotype but not phasing
    INCONSISTENT_PHASE=3
    #no call made (variantCall)
    UNCALLED=4
    #no call made (otherVariantCall)
    MISSING=5
    #calls a different allele
    DIFFERENT_ALLELE=6
    #suggests homozygous alternative (same non-reference alleles)
    POTENTIAL_HOMALT=7
    #suggests heterozygous (one reference allele)
    POTENTIAL_HET=8
    #suggests heterozygous (two different non-reference alleles)
    POTENTIAL_BIALLELIC=9
    
def agree_genotype(result):
    return result == Result.AGREE or \
           result == Result.AGREE_UNPHASED or \
           result == Result.INCONSISTENT_PHASE
def disagree_allele_count(result):
    return result == Result.POTENTIAL_HOMALT or \
           result == Result.POTENTIAL_HET or \
           result == Result.POTENTIAL_BIALLELIC
def disagree_genotype(result):
    return result == Result.DIFFERENT_ALLELE or \
           disagree_allele_count(result)
           
class Confidence(Enum):
    UNKNOWN=0
    CONSENSUS=1
    LIKELY_CANDIDATE=3
    CANDIDATE=2
    UNLIKELY_CANDIDATE=4
    SPURIOUS=5
    CALL=6

class Issues(Enum):
    INCOMPLETE_SUPPORT=1
    LOW_CONFIDENCE_CALL=2
    GENOTYPE_UNCERTAINTY=3
    QUESTIONABLE_PHASE=5

def compare(variantCall, otherCall):
    
    if variantCall is None and otherCall is None:
        return Result.AGREE
    if variantCall is None and otherCall is not None:
        return Result.UNCALLED
    if variantCall is not None and otherCall is None:
        return Result.MISSING

    # neither calls are None
    alleleA,alleleB = variantCall.alleleA,variantCall.alleleB
    otherAlleleA,otherAlleleB = otherCall.alleleA,otherCall.alleleB

    possibleAlleles = []
    if alleleA is not None and alleleA.allele not in possibleAlleles:
        possibleAlleles.append(alleleA.allele)
    if alleleB is not None and alleleB.allele not in possibleAlleles:
        possibleAlleles.append(alleleB.allele)
    if otherAlleleA is not None and otherAlleleA.allele not in possibleAlleles:
        possibleAlleles.append(otherAlleleA.allele)
    if otherAlleleB is not None and otherAlleleB.allele not in possibleAlleles:
        possibleAlleles.append(otherAlleleB.allele)

    def encode(call):
        a,b,p = 0, 0, 0
        if call.phased is not None and call.phased: p = 1
        A, B = call.alleleA, call.alleleB

        if A is not None:
            for i,allele in enumerate(possibleAlleles):
                if A.allele == allele:
                    a = i+1 ; break
        if B is not None:
            for i,allele in enumerate(possibleAlleles):
                if B.allele == allele:
                    b = i+1 ; break  
        
        return (a,b,p)

    v1, v2 = encode(variantCall), encode(otherCall)
    
    result = None    

    # agree homozygous
    if (v1[0],v1[1]) == (v2[0],v2[1]) and (v1[0] == v1[1]):
        result = Result.AGREE
    # agree but maybe unphased
    elif v1 == v2:
        result = Result.AGREE if v1[2] == 1 else  Result.AGREE_UNPHASED
    # same alleles
    elif (v1[0], v1[1]) == (v2[0], v2[1]) or \
        (v1[1], v1[0]) == (v2[0], v2[1]):

        # one or both unphased
         if v1[2] != 1 or v2[2] != 1:
             result = Result.AGREE_UNPHASED
    
         # both phased
         elif v1[2] == 1 and v2[2] == 1:
             result = Result.INCONSISTENT_PHASE

    # het disagreement
    elif v1[0] == 0 or v1[1] == 0:
        alt1 = v1[0] if v1[1] == 0 else v1[1]
        
        #alt does not match
        if v2[0] != alt1 and v2[1] != alt1:
             result = Result.DIFFERENT_ALLELE
        #alt matches both
        elif v2[0] == alt1 and v2[1] == alt1:
             result = Result.POTENTIAL_HOMALT
        #alt matches one (but v1 != v2)
        elif v2[0] == alt1 or v2[1] == alt1:
             result = Result.POTENTIAL_BIALLELIC

    # hom disagreement
    elif v1[0] == v1[1]:
        alt1 = v1[0]
        
        #alt does not match
        if v2[0] != alt1 and v2[1] != alt1:
             result = Result.DIFFERENT_ALLELE
        # v2 is het and alt matches one
        elif (v2[0] == 0 or v2[1] == 0):
             result = Result.POTENTIAL_HET
        # v2 is biallelic and matches one
        elif (v2[0] != 0 and v2[1] != 0) and (v2[0] == alt1 and v2[1] == alt1):
             result = Result.POTENTIAL_BIALLELIC 

    # v1 is biallelic
    else:
        result = Result.DIFFERENT_ALLELE
    
    return result
    
'''
class VariantCandidate:
    def __init__(self):
        self.primaryCall = None
        self.calls = dict()
        self.confidence = Confidence.UNKNOWN
        self.issues = []

    def add_call(self, name, variantCall):
        self.calls[name] = variantCall
        if self.primaryCall is None: self.primaryCall = variantCall
    def set_confidence(self, confidenceLevel):
        self.confidence=confidenceLevel
    def get_confidence(self):
        return self.confidence

    def add_issue(self, issue):
        self.issues.append(issue)
    def remove_issue(self, issue):
        if issue in self.issues: self.issues.remove(issue)
        
    def to_variant_call(self, ps=None):
        alleleA = self.primaryCall.get_alleleA(clean=True)
        alleleB = self.primaryCall.get_alleleB(clean=True)
        phased = self.primaryCall.is_phased()
                
        return VariantCall(alleleA, alleleB, phaseSet=ps, phased=phased)

    def __repr__(self):
        return self.get_confidence().name + ": " +  self.primaryCall.__repr__()
    def __str__(self):
       return "\n".join([call.__repr__() for call in self.calls])
'''

class Allele:
    def __init__(self, chrom, pos, ref, alt, vid="."): 
        self.chrom = chrom
        self.pos = pos
        self.ref = ref.upper()
        self.alt = alt.upper()
        self.vid = vid 
        
        self.althash = hash("-".join(str(x) for x in [chrom, pos, ref, alt]))
        self.refhash = hash("-".join(str(x) for x in [chrom, pos, ref]))

        self.homopolymer = None
        
    def get_chrom(self): return self.chrom
    def get_pos(self):   return self.pos
    def get_ref(self):   return self.ref
    def get_alt(self):   return self.alt

    def is_snp(self):
        return len(self.alt) == 1 and len(self.ref) == 1
    def is_insertion(self):
        return len(self.ref) == 1 and len(self.alt) > 1
    def is_deletion(self):
        return len(self.ref) > 1 and len(self.alt) == 1
    
    def shift(self, shift, newChrom=None):
        if newChrom is not None: self.chrom = newChrom
        self.pos += shift
    
    def __eq__(self, other): 
        return (self.chrom == other.chrom) and (self.pos == other.pos) and \
               (self.ref == other.ref) and (self.alt == other.alt)
                 
    def __len__(self): return 1
    
    def __hash__(self):
        return self.althash

    def get_vcf_cols(self):
        cols = [self.chrom, self.pos+1, self.vid, self.ref, self.alt]
        #note: vcf is 1-based
        return [str(x) for x in cols] 
     
    def __repr__(self):
        return "\t".join(self.get_vcf_cols())
    def __str__(self):
       return "\t".join(self.get_vcf_cols())

REF_ALLELE = Allele("__REF__", 0, "N", "N")

class VariantCall:
    
    def __init__(self, allele, tags=dict(), qual=".", info=".", filt="PASS"):
        self.pos = allele.get_pos()
        self.chrom = allele.get_chrom()
        self.ref = allele.get_ref()

        self.tags = tags
        self.qual = qual
        self.info = info
        self.filter = filt

        self.alleles = [REF_ALLELE, allele]
        self.call = [REF_ALLELE, allele]
        self.phaseSet = "."

    def set_call(self, alleleA=None, alleleB=None, ps="."):
        if alleleA is None: alleleA = REF_ALLELE
        if alleleB is None: alleleB = REF_ALLELE
        
        if alleleA not in self.alleles: self.alleles.append(alleleA)
        if alleleB not in self.alleles: self.alleles.append(alleleB)

        self.call = [alleleA, alleleB]

        self.phaseSet = ps
        
    def add_allele(self, allele, asCall=True):
        if allele == REF_ALLELE: return
        
        if allele.get_pos() != self.pos: 
            print("position issue while adding allele:\n",
                  self.to_vcf_line(), "\n", allele)
        if allele.get_chrom() != self.chrom: 
            print("chromosome issue while adding allele:\n",
                  self.to_vcf_line(), "\n", allele)
        if allele.get_ref() != self.ref: 
            print("ref issue while adding allele:\n",
                  self.to_vcf_line(), "\n", allele)
        
        if allele not in self.alleles:
            self.alleles.append(allele)
        
        if asCall:
            if self.call[0] == REF_ALLELE:
               self.call = [allele, self.call[1]] 
            if self.call[1] == REF_ALLELE:
               self.call = [self.call[0], allele] 
        
    def add_alleles(self, otherVariantCall, asCall=True):
        for allele in otherVariantCall.alleles:
            self.add_allele(allele, asCall)

    def flip_alleles(self):
        self.call = [self.call[0], self.call[1]]
    
    def get_pos(self): return self.pos
    def get_chrom(self): return self.chrom
    
    def is_heterozygous(self):
        return self.call[0] != self.call[1]    
    def is_homozygous_alt(self):
        return self.call[0] == self.call[1] and self.call[0] != REF_ALLELE    
    def is_homozygous_ref(self):
        return self.call[0] == REF_ALLELE and self.call[1] == REF_ALLELE
    def is_biallelic_alt(self):
        return self.call[0] != self.call[1] and self.call[0] != REF_ALLELE and self.call[1] != REF_ALLELE
       
    def is_snp(self):
        if self.is_homozygous_ref(): return False
        for call in self.call:
            if call != REF_ALLELE and not call.is_snp(): return False
        return True
    
    def is_insertion(self):
        if self.is_homozygous_ref(): return False
        for call in self.call:
            if call != REF_ALLELE and not call.is_insertion(): return False
        return True

    def is_deletion(self):
        if self.is_homozygous_ref(): return False
        for call in self.call:
            if call != REF_ALLELE and not call.is_deletion(): return False
        return True

    def main_allele(self):
        return self.alleles[1]
    def update_main_allele(self, newAllele):
        return self.alleles[1]

    def is_phased(self):
        return self.phaseSet != "."
    def get_phase_set(self):
        return self.phaseSet
    def set_phased(self, ps="0"):
        self.phaseSet = ps

    def get_gt_string(self):
        gt = []
        
        for gtAllele in self.call:
            for i,allele in enumerate(self.alleles):
                if gtAllele == allele:
                    gt.append(str(i))
        
        sep = "/" if self.phaseSet == "." else "|"
        
        return sep.join(gt)
    
    def get_format_string(self):
        frmt = ["GT", "PS"]
        data = [self.get_gt_string(), str(self.phaseSet)]
        
        for tag in self.tags:
            if tag not in frmt:
                frmt.append(tag)
                data.append(self.tags[tag])
        
        sep = ":"
        return (sep.join(frmt), sep.join(data))
    
    def to_vcf_line(self):
        #if self.is_homozygous_ref():
        #    return ""

        alts = [allele.get_alt() for allele in self.alleles if allele != REF_ALLELE]
        alt = ",".join(alts)     
        vids = [allele.vid for allele in self.alleles if allele != REF_ALLELE]
        vid = ",".join(list(set(vids)))
        
        #note: vcf is 1-based
        vcfCols = [self.chrom, self.pos+1, vid, self.ref, alt]

        frmt, data = self.get_format_string()

        vcfCols = vcfCols + [self.qual, self.filter, self.info, frmt, data]
        vcfLine = "\t".join([str(x) for x in vcfCols])   
        
        return vcfLine
        
    def shift(self, shift, newChrom=None):
        self.pos += shift
        if newChrom is not None:
            self.chrom = newChrom

        for allele in self.alleles:
            if allele != REF_ALLELE: 
                allele.shift(shift, newChrom)
    
    def __repr__(self):
        return self.to_vcf_line()
    def __str__(self):
       return self.to_vcf_line()


class VariantCallSet:
    def __init__(self):
        self.calls = dict()
        self.filtered = []
        self.defaultChrom = None
        
    def add(self, variantCall):
        chrom, pos = variantCall.get_chrom(), variantCall.get_pos() 
        if self.defaultChrom is None: self.defaultChrom = chrom

        if chrom in self.calls:
            if pos in self.calls[chrom]:
                self.calls[chrom][pos].add_alleles(variantCall)
            else:
                self.calls[chrom][pos] = variantCall
        else:
            
            self.calls[chrom] = dict()
            self.calls[chrom][pos] = variantCall
            
    def add_all(self, otherCallset):
        for chrom in otherCallset.calls:
            for pos in otherCallset.calls[chrom]:
                self.add(otherCallset.get(pos, chrom))

    def get(self, pos, chrom=None):
        if chrom is None: chrom = self.defaultChrom
        
        if chrom not in self.calls: return None
        if pos not in self.calls[chrom]: return None        
        return self.calls[chrom][pos]
    
    def drop_outside(self, region):
        self.defaultChrom = region.chrom
        for chrom in self.calls:
            if chrom != region.chrom:
                self.calls.pop(chrom)
        
        for pos in self.get_positions(region.chrom):
            if not region.contains(region.chrom, pos):
                self.calls.pop(pos)
        
    def pop(self, pos, chrom=None):
        if chrom is None: chrom = self.defaultChrom
        
        if chrom not in self.calls: return None
        if pos not in self.calls[chrom]: return None
        return self.calls[chrom].pop(pos)

        if pos in self.calls:
            return self.calls.pop(pos)
        else: return None
    
    def remove(self, variantCall):
        self.pop(variantCall.get_pos(), variantCall.get_chrom())

    def get_positions(self, chrom=None):
        if chrom is None: chrom = self.defaultChrom
        return sorted(self.calls[chrom])
    
    def get_chroms(self):
        return sorted(self.calls)

    def get_variant_calls(self, chrom=None, allChrom=False):
        if chrom is None: chrom = self.defaultChrom
        callList = []

        for ch in self.get_chroms():
            if not allChrom and ch != chrom: continue
            callList.extend([self.calls[ch][pos] for pos in self.get_positions(ch)])
        return callList
    
    def get_range(self, chrom=None):
        if chrom is None: chrom = self.defaultChrom
        positions = self.get_positions(chrom)
        return (positions[0], positions[-1])

    def flip_alleles(self, chrom=None, allChrom=False):
        if chrom is None: chrom = self.defaultChrom
        
        for ch in self.calls:
            if not allChrom and ch != chrom: continue
            for pos in self.calls[ch]:
                self.calls[ch][pos].flip_alleles()            
            
    def set_all_phased(self, chrom=None, allChrom=False, ps="0"):
        if chrom is None: chrom = self.defaultChrom
        
        for ch in self.calls:
            if not allChrom and ch != chrom: continue
            for pos in self.calls[ch]:
                self.calls[ch][pos].set_phased(ps)       
    
    def shift_all(self, shift, newChrom=None, chrom=None):
        if chrom is None: chrom = self.defaultChrom

        shiftedDict = dict()
        for pos in self.get_positions(chrom):
            vc = self.calls[chrom].pop(pos)
            vc.shift(shift, newChrom)
            shiftedDict[vc.get_pos()] = vc
        self.calls[chrom] = shiftedDict
    
    def write_vcf(self, vcfName, compressed=False):
        tab, nl = '\t', '\n'

        vcf = open(vcfName, "w+")

        vcf.write("##fileformat=VCFv4.2" + nl)
        
        header = tab.join(["#CHROM", "POS", "ID", "REF", "ALT",
                           "QUAL",  "FILTER", "INFO", "FORMAT", "SAMP"])
    
        vcf.write(header + nl)
    
        for variantCall in self.get_variant_calls(allChrom=True):
            vcf.write(variantCall.to_vcf_line() + nl)
                
        vcf.close()
                        
        if compressed:
            subprocess.call(["bgzip", "-f", vcfName])
            subprocess.call(["tabix", vcfName + ".gz"])
            return vcfName + ".gz"       
            #todo: delete unzipped????
        
        return vcfName
    
    def __len__(self): return sum([len(self.calls[ch]) for ch in self.calls])

    def __repr__(self):
        return "Variant set with " + str(len(self)) + " variants"
    
    def __str__(self):
         return "\n".join([v.__repr__() for v in self.get_variant_calls(allChrom=True)])

def split_by_ps(variantSet):
    calls = variantSet.get_variant_calls()
    splitSets = dict()

    for call in calls:
        ps = call.get_phase_set()
        if ps not in splitSets:
            splitSets[ps] = VariantCallSet()            
        splitSets[ps].add(call)

    return splitSets

def vcf_to_variantcallset(vcfName, region=None, tigId=None, liftoverRegion=None, ignoreNonVariant=True):
    variantSet = VariantCallSet()
    with open(vcfName) as vcf:
        line = True
        while line:
            line = vcf.readline()
            if line is None or len(line) < 1 or line[0] == '#': continue
            cols = (line.rstrip("\n\r")).split('\t')
           
            chrom,pos,ref,alt = cols[0],int(cols[1]),cols[3],cols[4]

            if tigId is not None: 
                if chrom != tigId: continue
            if region is not None:
                if pos < region.start or pos > region.end: continue
                if chrom != region.chrom: continue
            
            if ignoreNonVariant and ref.upper() == alt.upper():
                continue
            
            shift=0
            if liftoverRegion is not None:
                chrom = liftoverRegion.chrom
                shift = liftoverRegion.start 
            
            alts = alt.split(",")
            alleles = [REF_ALLELE]

            for alt in alts:
                #note:vcf are 1-based!
                allele = Allele(chrom, pos-1 + shift, ref, alt)
                alleles.append(allele)
            
            qual,filt,info = cols[5],cols[6],cols[7]
            tags = dict()
            gtA,gtB = None, None
            ps = "."
            
            if len(cols) > 8:
                frmt, data = cols[8], cols[9]
                
                for f,s in zip(frmt.split(':'), data.split(':')):
                    
                    tags[f] = s
                    
                    if f == "GT":
                        gtA = None if s[0] == "." else int(s[0])
                        gtB = None if s[-1] == "." else int(s[-1])
                        #phased = ('|' in s)
                    elif f == "PS":
                        ps = s
            
            if gtA is None and gtB is None: 
                gtA=0 ; gtB=1

            variantCall = VariantCall(alleles[1], tags=tags, qual=qual, info=info, filt=filt)
            variantCall.set_call(alleleA=alleles[gtA], alleleB=alleles[gtB], ps=ps)
            variantSet.add(variantCall)
            
    vcf.close()    
    return variantSet

def collapse_variants_between(variantSetA, variantSetB, seq, collapseRange=2):
    
    #B comes before A
    def collapse(A, posA, B, posB, dist):
        
        a = A.main_allele()
        b = B.main_allele()

        anchorEnd = max(posB + len(b.ref), posA + len(b.ref))
        anchor = seq[posB:anchorEnd]
        newA = Allele(a.chrom, posB, 
                       anchor[:dist] + a.ref + anchor[dist + len(a.ref):], 
                       anchor[:dist] + a.alt + anchor[dist + len(a.ref):])
        newB = Allele(b.chrom, posB, 
                       b.ref + anchor[len(b.ref):], 
                       b.alt + anchor[len(b.ref):])
        
        
        return (newA, newB)
    
    for chrom in set(variantSetA.get_chroms() + variantSetB.get_chroms()):
        done = dict()
        complete=False
        while(not complete):
            
            positions = variantSetA.get_positions(chrom)
            
            for posA in positions:
    
                if posA in done: continue
                flag, complete = False, False
                
                if variantSetB.get(posA, chrom) is not None and \
                    variantSetA.get(posA, chrom).ref != variantSetB.get(posA, chrom).ref:
                    A, B = variantSetA.pop(posA, chrom), variantSetB.pop(posA, chrom)
                    newA, newB = collapse(A, posA, B, posA, 0)
                    variantSetA.add(newA)
                    variantSetB.add(newB)
                    break
                
                for i in range(collapseRange):
                    posB = posA-(i+1)
                    if variantSetB.get(posB, chrom) is not None:
                        flag = True
                        A, B = variantSetA.pop(posA, chrom), variantSetB.pop(posB, chrom)
                        newA, newB = collapse(A, posA, B, posB, i+1)
                        variantSetA.add(newA)
                        variantSetB.add(newB)
                if flag: break
    
                for i in range(collapseRange):
                    posB = posA+(i+1)
                    if variantSetB.get(posB, chrom) is not None:
                        flag = True
                        A, B = variantSetA.pop(posA, chrom), variantSetB.pop(posB, chrom)
                        newB, newA = collapse(B, posB, A, posA, i+1)
                        variantSetA.add(newA)
                        variantSetB.add(newB)                   
                if flag: break
            
            
            
                done[posA] = True
                complete = True
    
    return(variantSetA, variantSetB)


def combine_as_genotypes(variantSetA, variantSetB, phased=False):
    variantSetC = copy.deepcopy(variantSetA)
    
    variantSetB.flip_alleles()
    variantSetC.add_all(variantSetB)
    variantSetC.set_all_phased()
    
    return variantSetC


def create_variant_call(chrom, pos, ref, alt):
    allele = Allele(chrom, pos, ref, alt)
    variantCall = VariantCall(allele)
    return variantCall   
    
    
    
'''
class Variant:
    def __init__(self, chrom, pos, ref, alt, annotation=""): 
        
        self.chrom = chrom
        self.pos = pos
        self.ref = ref.upper()
        self.alt = alt.upper()
        self.otherAlt = []

        self.annotation = annotation

        self.vid = "."
        self.qual= "."
        self.info = "."
        self.filter = "PASS"
        
        self.tags = dict()
        self.gt_a = None
        self.gt_b = None
        self.ps = None
        self.homopolymer = None

    def set_gt(self, a, b, ps="."):
        self.gt_a = a
        self.gt_b = b
        self.ps = ps

    def get_allele(self, gt):
        if gt == 0: return self.ref
        if gt == 1: return self.alt
        if gt > 1: return self.otherAlt[gt-2]
        return None
    
    def get_gt(self, phased=False):
        if self.gt_a is None or self.gt_b is None: 
            return None
        
        if phased:
            return (self.gt_a, self.gt_b)
        
        return self.gt_a + self.gt_b

    def add_alt(self, newAlt):
        self.otherAlt.append(newAlt)

    def get_gt_string(self):
        if self.gt_a is None and self.gt_b is None:
            return "0/1"
            
        if self.gt_a is not None and self.gt_b is not None:
            phased = "/"
            if self.ps is not None and self.ps != ".": phased = "|"
            
            return str(self.gt_a) + phased + str(self.gt_b)
        return None
    
    def get_format_string(self):
        frmt="GT"
        data=self.get_gt_string()
        if self.ps is not None:
            frmt = frmt + ":PS"
            data = data + ":" + str(self.ps)

        if self.homopolymer is not None:
            frmt = frmt + ":" + HOMOPOLYMER_TAG
            data = data + ":" + self.homopolymer

        return (frmt, data)

    
    def to_vcf_line(self):
        frmt, data = self.get_format_string()
        
        alt = self.alt + ("," if len(self.otherAlt) > 0 else "") + \
            ",".join([a for a in self.otherAlt])
        
        #note: vcf is 1-based
        vcfList = [self.chrom, self.pos+1, self.vid, self.ref, alt, 
                   self.qual, self.filter, self.info, frmt, data]
        
        vcfLine = "\t".join([str(x) for x in vcfList])   
        return vcfLine
    
    def is_snp(self):
        if not len(self.ref) == 1: return False
        if not len(self.alt) == 1: return False
        for alt in self.otherAlt: 
            if not len(alt) == 1: return False
        return True
    
    def is_insertion(self):
        #multiple alt???
        return len(self.ref) == 1 and len(self.alt) > 1

    def is_deletion(self):
        #multiple alt???
        return len(self.ref) > 1 and len(self.alt) == 1
    
    def variant_compare(self, other):
        frmt="GT"
        data=self.get_gt_string()
        if self.ps is not None:
            frmt = frmt + ":PS"
            data = data + ":" + str(self.ps)

        if self.homopolymer is not None:
            frmt = frmt + ":" + HOMOPOLYMER_TAG
            data = data + ":" + self.homopolymer

        return (frmt, data)

    
    def get_copy(self):
        return copy.deepcopy(self)
    
    def __eq__(self, other): 
        return (self.chrom == other.chrom) and (self.pos == other.pos) and \
               (self.ref == other.ref) and (self.alt == other.alt)
                 
    def __len__(self): return 1
     
    def __repr__(self):
        return "\t".join([str(x) for x in (self.chrom, self.pos, self.ref, self.alt)])
    def __str__(self):
        return self.to_vcf_line()
    
    

class VariantSet:
    def __init__(self): 
        self.vs = dict()
        self.filtered = []

    def get(self, pos):
        if pos in self.vs:
            return self.vs[pos]
        else:
            return None
        
    def add(self, variant):
        if variant.pos in self.vs:
            print("Warning: variant position already exists in variant set:\n" + str(variant))
        self.vs[variant.pos] = variant
        
    def remove(self, variant):
        if variant.pos in self.vs:
            if self.vs[variant.pos] == variant:
                self.vs.pop(variant.pos)
        
    def pop(self, pos):
        if pos in self.vs:
            return self.vs.pop(pos)
        else: return None

    def get_positions(self):
        return sorted(self.vs)
    def get_variants(self):
        return [self.vs[pos] for pos in self.get_positions()]
       
    def get_range(self):
        positions = self.get_positions()
        return (positions[0], positions[-1])

    def filter_by_tag(self, tag, value, eq=False, lt=False, gt=False, keepMissing=False):
        if not eq and not lt and not gt:
            return
        
        for pos in self.get_positions():
            variant = self.vs[pos]
            if not tag in variant.tags or variant.tags[tag] == ".":
                if not keepMissing: self.filtered.append(self.pop(pos))
                continue
                
            if eq and variant.tags[tag] == value:
                continue
                
            if lt and int(variant.tags[tag]) < value:
                continue
            if gt and int(variant.tags[tag]) > value:
                continue
            
            self.filtered.append(self.pop(pos))
    
    def filter_homopolymer(self):
        self.filter_by_tag(HOMOPOLYMER_TAG, HOMOPOLYMER_NO, eq=True)
    
    def filter_nonhet(self):
        for pos in self.get_positions():
            variant = self.vs[pos]
            gt = variant.get_gt()
            if gt is None or gt != 1:
                self.filtered.append(self.pop(pos))    

    def filter_nonsnp(self):
        for pos in self.get_positions():
            variant = self.vs[pos]
            if not variant.is_snp():
                self.filtered.append(self.pop(pos))    

    def unfilter_snps(self):
        newFiltered = []
        
        while len(self.filtered) > 0:
            variant = self.filtered.pop()
            if variant.is_snp():
                self.vs[variant.pos] = variant
            else: newFiltered.append(variant)
        self.filtered = newFiltered
    
    def unfilter_matching_pos(self, otherSet):
        newFiltered = []
        
        while len(self.filtered) > 0:
            variant = self.filtered.pop()
            if variant.pos in otherSet.vs:
                self.vs[variant.pos] = variant
            else: newFiltered.append(variant)
        self.filtered = newFiltered
    
    def write_vcf(self, vcfName, hetsOnly=False, compressed=False):
    
        tab, nl = '\t', '\n'
        header = tab.join(["#CHROM", "POS", "ID", "REF", "ALT",
                           "QUAL",  "FILTER", "INFO", "FORMAT", "SAMP"])
        
        vcf = open(vcfName, "w+")
        vcf.write("##fileformat=VCF" + nl)
        vcf.write(header + nl)
        
        for variant in self.get_variants():
            if not hetsOnly or variant.get_gt() == 1:
                vcf.write(variant.to_vcf_line() + nl)
                
        vcf.close()
                        
        if compressed:
            bgzipped = tools.bgzip(vcfName)
            #todo: delete unzipped????
            return bgzipped
        
        return vcfName

    def summarize(self, printLines=True):
        
        nvariants = len(self.vs)
        if nvariants == 0:
            nope = "No variants."
            if printLines: 
                print(nope)
                return nope
            return [nope]
    
        nphased = sum([v.ps is not None and v.ps != "." for v in self.get_variants()])
        homref = sum([v.gt_a == 0 and v.gt_b == 0 for v in self.get_variants()])        
        hetA = sum([v.gt_a > 0 and v.gt_b == 0 for v in self.get_variants()]) 
        hetB = sum([v.gt_a == 0 and v.gt_b > 0 for v in self.get_variants()]) 
        homalt = sum([v.gt_a > 0 and v.gt_b > 0 for v in self.get_variants()]) 
        
        longestRef = None
        longestAlt = None
        
        for v in self.get_variants():
            if longestRef is None or v.ref > longestRef.ref:
                longestRef = v
            if longestAlt is None or v.alt > longestAlt.alt:
                longestAlt = v

        def percent(x):
            return str(round(100*x/nvariants,2)) + "%"

        lines = [str(nvariants) + " variants, " + str(nphased) + " phased.",
                 percent(homref) + " 0/0\t" + percent(hetA + hetB) + " 0/1\t " + percent(homalt) + " 1/1",
                 percent(hetA) + " on haplotype A\t" + percent(hetB) + " on haplotype B\t",
                 "longest ref allele:",
                 longestRef.to_vcf_line(),
                 "longest alt allele:",
                 longestAlt.to_vcf_line()]

        summary = '\n'.join(lines)

        if printLines: 
            print(summary)
            return summary
        return lines
    
    def __len__(self): return len(self.vs)

    def __repr__(self):
        return "Variant set with " + str(len(self.vs)) + " variants"
    
    def __str__(self):
        return "\n".join([v.__repr__() for pos,v in sorted(self.vs.items())])



def combine_as_genotypes(variantSetA, variantSetB, seq, collapseRange=2):

    variantSetA = collapse_variants_self(variantSetA, seq, collapseRange)
    variantSetB = collapse_variants_self(variantSetB, seq, collapseRange)
    variantSetA, variantSetB = collapse_variants_between(variantSetA, variantSetB, seq, collapseRange=2)
    
    variantSetGT = VariantSet()
    
    for pos in variantSetA.get_positions():
        A = variantSetA.pop(pos)
        B = variantSetB.pop(pos)
                
        if B is not None:            
            if A == B:
                gtB = 1
            elif A.ref == B.ref:
                gtB = 2
                A.add_alt(B.alt)
            else:
                print("Warning: variants start at same position but have different REF:\n" + \
                      str(A) + "\n" + str(B))
        else:
            gtB = 0
            
        A.set_gt(1, gtB, ps="0")
        variantSetGT.add(A) 
        
    for pos in variantSetB.get_positions():
        B = variantSetB.pop(pos)
        B.set_gt(0, 1, ps="0")
        variantSetGT.add(B) 

    return variantSetGT
        
    
def collapse_variants_self(variantSet, seq, collapseRange=2):
    
    done = dict()
    complete=False
    while(not complete):
        
        positions = variantSet.get_positions()
        
        for pos1 in positions:
            if pos1 in done: continue
            flag, complete = False, False

            for i in range(collapseRange):
                pos2 = pos1-(i+1)
                if variantSet.get(pos2) is not None:
                    flag = True
                    V1 = variantSet.pop(pos1)
                    V2 = variantSet.pop(pos2)
                    anchor = seq[pos2:pos1]
                    ref = V2.ref + anchor[len(V2.ref):] + V1.ref
                    alt = V2.alt + anchor[len(V2.ref):] + V1.alt
                    
                    V = Variant(V1.chrom, pos2, ref, alt)
                    variantSet.add(V)
            if flag: break
        
            for i in range(collapseRange):
                pos2 = pos1+(i+1)
                if variantSet.get(pos2) is not None:
                    flag = True
                    V1 = variantSet.pop(pos1)
                    V2 = variantSet.pop(pos2)
                    anchor = seq[pos1:pos2]
                    ref = V1.ref + anchor[len(V1.ref):] + V2.ref
                    alt = V1.alt + anchor[len(V1.ref):] + V2.alt

                    V = Variant(V1.chrom, pos1, ref, alt)
                    variantSet.add(V)
            if flag: break
        
            done[pos1] = True
            complete = True
    
    return(variantSet)


def collapse_variants_between(variantSetA, variantSetB, seq, collapseRange=2):
    
    #B comes before A
    def collapse(A, posA, B, posB, dist):
                
        anchorEnd = max(posB + len(B.ref), posA + len(A.ref))
        anchor = seq[posB:anchorEnd]
        newA = Variant(A.chrom, posB, 
                       anchor[:dist] + A.ref + anchor[dist + len(A.ref):], 
                       anchor[:dist] + A.alt + anchor[dist + len(A.ref):])
        newB = Variant(B.chrom, posB, 
                       B.ref + anchor[len(B.ref):], 
                       B.alt + anchor[len(B.ref):])
        return (newA, newB)
    
    done = dict()
    complete=False
    while(not complete):
        
        positions = variantSetA.get_positions()
        
        for posA in positions:

            if posA in done: continue
            flag, complete = False, False
            
            if variantSetB.get(posA) is not None and \
                variantSetA.get(posA).ref != variantSetB.get(posA).ref:
                A, B = variantSetA.pop(posA), variantSetB.pop(posA)
                newA, newB = collapse(A, posA, B, posA, 0)
                variantSetA.add(newA)
                variantSetB.add(newB)
                break
            
            for i in range(collapseRange):
                posB = posA-(i+1)
                if variantSetB.get(posB) is not None:
                    flag = True
                    A, B = variantSetA.pop(posA), variantSetB.pop(posB)
                    newA, newB = collapse(A, posA, B, posB, i+1)
                    variantSetA.add(newA)
                    variantSetB.add(newB)
            if flag: break

            for i in range(collapseRange):
                posB = posA+(i+1)
                if variantSetB.get(posB) is not None:
                    flag = True
                    A, B = variantSetA.pop(posA), variantSetB.pop(posB)
                    newB, newA = collapse(B, posB, A, posA, i+1)
                    variantSetA.add(newA)
                    variantSetB.add(newB)                   
            if flag: break
        
        
        
            done[posA] = True
            complete = True
    
    return(variantSetA, variantSetB)

def set_homopolymer_status(variantSet, seqData, minRepeat=1):
    
    for variant in variantSet.get_variants():
        bp = None
        unique = False
                       
        if variant.is_deletion():
            bp = variant.ref[1]
            unique = len(set(variant.ref[1:])) == 1
        if variant.is_insertion():
            bp = variant.alt[1]
            unique = len(set(variant.alt[1:])) == 1
        
        if bp is None or unique == False:
            variant.homopolymer=HOMOPOLYMER_NO
            continue
        
        #for debugging
        print("------")
        print(variant.ref, variant.alt)
        print(bp)
        print(seqData[variant.chrom][variant.pos: variant.pos + 10])
        start = variant.pos+1
        print(seqData[variant.chrom][start:start + minRepeat])
        
        start = variant.pos+1
        refBps = seqData[variant.chrom][start:start + minRepeat]
        
        flag = False
        for refBp in refBps:
            if bp != refBp: 
                variant.homopolymer=HOMOPOLYMER_NO
                flag = True
                break
        if flag: continue

        variant.homopolymer=HOMOPOLYMER_YES
        
        #for debugging
        #print("YES")

    for variant in variantSet.get_variants():
        variant.tags[HOMOPOLYMER_TAG] = str(variant.homopolymer)
        

'''
    