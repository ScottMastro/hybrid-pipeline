import pandas as pd

'''    
lrdict = scaffolder.linkedReadsDict(rids, param.arks_output)
#lrdf = scaffolder.linkedReadsDF(rids, lengthData, param.arks_output)
unitigs = scaffolder.unitigsDict(param.unitigsBed)
'''

def checkLR(p1, lrdict):
#checking how often path-switches are found in linked reads (when canu tig changes)
    distrib = dict()
    distrib[0] = 0
    dists = dict()
    dists[0] = list()
    LRtigs = []
    pathTigs = []
    notPathTigs = []
    for p in p1:
        #if path_helper.path_length(p) < 10000000:
            #continue
        for i in range(0, len(p)-1):
            if p[i].rid not in pathTigs:
                pathTigs.append(p[i].rid)
            if p[i+1].rid not in pathTigs:
                pathTigs.append(p[i+1].rid)
            if p[i].rid != p[i+1].rid and p[i].rid != 'NNN' and p[i+1].rid != 'NNN':
                tig1 = p[i].rid
                tig2 = p[i+1].rid
                pos1 = p[i].after_pos()
                
                pos2 = p[i+1].before_pos()
                dist = abs(pos2 - pos1)
                if tig1 in lrdict:
                    a = lrdict[tig1]
                    if tig2 in a:
                        n = lrdict[tig1][tig2]
                        if n not in distrib:
                            distrib[n] = 1
                            dists[n] = list()
                            dists[n].append(dist)
                        else:
                            distrib[n] += 1
                            dists[n].append(dist)
                    else:
                        distrib[0] += 1
                        dists[0].append(dist)
                        #print('tig2: '+tig2)
                        #print(p[i])
                        #print(p[i+1])
                        #print('===============')
                else:
                    distrib[0] += 1
                    dists[0].append(dist)
                    #print('tig1: '+tig1)

                
    distribLR = dict()
    for tig in lrdict:
        if tig not in LRtigs:
            LRtigs.append(tig)
        for tig2 in lrdict[tig]:
            if lrdict[tig][tig2] not in distribLR:
                distribLR[lrdict[tig][tig2]] = 1
            else:
                distribLR[lrdict[tig][tig2]] += 1
    for tig in LRtigs:
        if tig not in pathTigs:
            notPathTigs.append(tig)

def linkedReadsDF(canutigs, lengthData, arks_output):
    #Make df with linked reads, currently unused
    lrdf = pd.DataFrame(columns=['tig1', 'tig1len', 'dir2', 'tig2', 'tig2len', 'dir2', 'n_barcodes'])
    file = open(arks_output, 'r')
    i = 0
    freq_dict = dict()
    for tig in canutigs:
        freq_dict[tig] = 0 
    
    for line in file:
        sline = line.split('\t')
        first = sline[1]
        second = sline[2]
        tig1_r = first[1:] #tig1_r = tig1_renamed --> contig in the arks renamed fasta file
        tig2_r = second[1:]
        
        tig1 = canutigs[(int(tig1_r))-1] #finding the original contig from the renamed one
        tig1len = lengthData[tig1]
        dir1 = first[0]
        tig2 = canutigs[(int(tig2_r))-1]
        tig2len = lengthData[tig2]
        dir2 = second[0]
        bars = int(sline[3])
        
        freq_dict[tig1] += bars
        freq_dict[tig2] += bars
        
        lrdf.loc[i] = tig1, tig1len, dir1, tig2, tig2len, dir2, bars
        i += 1
    #tig1_freq and tig2_freq = total number of barcodes for this contig that connect to contigs other than the one in the current row
    list1 = []
    list2 = []
    for index, row in lrdf.iterrows():
        list1.append(freq_dict[row[0]]-row[6]) 
        list2.append(freq_dict[row[3]]-row[6])
    lrdf.insert(2, 'tig1_freq', list1)
    lrdf.insert(6, 'tig2_freq', list2)
    
    file.close()
    lrdf.to_csv('linked_reads.csv')
    
    return lrdf

def checkLinkedReads(path, nextPath, linkedReads, threshold):   
    tig1a = path[0].rid
    tig1b = path[-1].rid
    tig2a = nextPath[0].rid
    tig2b = nextPath[-1].rid
    
    if tig1a in linkedReads.keys():
        if tig2a in linkedReads[tig1a].keys():
            if linkedReads[tig1a][tig2a] >= threshold:
                return True
        if tig2b in linkedReads[tig1a].keys():
            if linkedReads[tig1a][tig2b] >= threshold:
                return True
    if tig1b in linkedReads.keys():
        if tig2a in linkedReads[tig1b].keys():
            if linkedReads[tig1b][tig2a] >= threshold:
                return True
        if tig2b in linkedReads[tig1b].keys():
            if linkedReads[tig1b][tig2b] >= threshold:
                return True
    return False

def unitigsDict(unitigs):
    utgdict_dict = dict()
    file = open(unitigs, 'r')
    for line in file:
        sline = line.split('\t')
        ctg_raw = sline[0]
        ctg = "tig"+ctg_raw[3:]+"_pilon_pilon"
        start = sline[1]
        end = sline[2]
        rng = (start, end)
        utg = sline[3]
        
        if ctg not in utgdict_dict.keys():
            utgdict = dict()
            utgdict[utg] = rng
            utgdict_dict[ctg] = utgdict
        else:
            utgdict_dict[ctg][utg] = rng
    return utgdict_dict

def linkedReadsDict(canutigs, arks_output):
    lrdict_dict = dict()
    file = open(arks_output, 'r')
    for line in file:
        sline = line.split('\t')
        first = sline[1]
        second = sline[2]
        tig1_r = first[1:] #tig1_r = tig1_renamed --> contig in the arks renamed fasta file
        tig2_r = second[1:]
        tig1 = canutigs[(int(tig1_r))-1] #finding the original contig from the renamed one
        tig2 = canutigs[(int(tig2_r))-1]
        bars = int(sline[3])
        
        if tig1 not in lrdict_dict.keys():
            lrdict = dict()
            lrdict[tig2] = bars
            lrdict_dict[tig1] = lrdict
        else:
            lrdict_dict[tig1][tig2] = bars
    return lrdict_dict