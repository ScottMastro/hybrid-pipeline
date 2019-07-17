from parameters import Parameters
import file_handler as io
import pandas as pd

def region_extract(seqData, contig, pos1, pos2):
    seq = seqData[str(contig)]
    region = seq[pos1-1:pos2-1]
    return region

def findOverlaps(paths_c, lengthData):
    from path_helper import path_overlap2

    file = open("overlaps.txt", "w")
    overlaps = {}
    for i in range(0, len(paths_c)):
        for j in range(0, len(paths_c)):
            if i >= j:
                continue
            overlap = path_overlap2(paths_c[i], paths_c[j], lengthData, True)
            #print(str(overlap) + ', ' + str(i) + ', ' + str(j)) 
            if overlap > 0:
                overlaps[str(i)+','+str(j)] = overlap
                print(str(i) + ', ' + str(j))
                file.write(str(i) + "\t" + str(j) + "\n")
    
    
    df = pd.DataFrame(columns=['contigs', 'pathID1', 'whichEnd1', 'posA', 'nforks1', 'pathID2', 'whichEnd2', 'posB', 'nforks2', 'repeat_size'])
    file = open("overlaps.txt", "r")
    i=0
    for line in file:
        splitline = line.split()
        path1 = int(splitline[0])
        nforks1 = len(paths_c[path1])
        path2 = int(splitline[1])
        nforks2 = len(paths_c[path2])
        whichEnd1, whichEnd2, tigs, posA, posB, sizeA, sizeB = path_overlap2(paths_c[path1], paths_c[path2], lengthData, False)
        
        for j in range(0, len(tigs)):
            dfwhichEnd1 = whichEnd1[j]
            dfwhichEnd2 = whichEnd2[j]
            dftigs = tigs[j]
            dfposA = posA[j]
            dfposB = posB[j]
            dfsizeA = sizeA[j]
            dfsizeB = sizeB[j]
            repeat = abs(dfsizeB - dfsizeA)
            df.loc[i] = dftigs, path1, dfwhichEnd1, dfposA, nforks1, path2, dfwhichEnd2, dfposB, nforks2, repeat
            i += 1
    file.close()
    df.to_csv('overlap_starts_ends.csv')
    
    return df

def findOverlaps2(df, lengthData, seqData):
#df2 = pd.DataFrame(columns=['region1', 'alignment1', 'region2', 'alignment2', 'region3', 'alignment3', 'region4', 'alignment4', 'region5', 'alignment5'])
    regionList = []
    for row in df.iterrows():
        tig = row[1]['contigs']
        if tig == "NNN":
            continue
        whichEnd1 = row[1]['whichEnd1']
        posA = row[1]['posA']
        posA1 = int(posA.split('-')[0])
        posA2 = int(posA.split('-')[1])
        whichEnd2 = row[1]['whichEnd2']
        posB = row[1]['posB']
        posB1 = int(posB.split('-')[0])
        posB2 = int(posB.split('-')[1])
        posList = [posA1, posA2, posB1, posB2]
        posList.sort()
        posList.append(max((posList[0]-50000), 1))
        posList.append(min((posList[3]+50000),lengthData[str(tig)]-1))
        posList.sort()
        
        
        for k in range(0, 5):
            dictkey = tig + ',' + str(posList[k]) + '-' + str(posList[k+1])
            dicttuple = dictkey, region_extract(seqData, tig, posList[k], posList[k+1])
            regionList.append(dicttuple)
            
        smallest = min(posA1, posA2, posB1, posB2)
        largest = max(posA1, posA2, posB1, posB2)
        if whichEnd1 == 'Head':
            if whichEnd2 == 'Head':
                
                print('Head Head')
            elif whichEnd2 == 'Middle':
                print('Head Middle')
            elif whichEnd2 == 'Tail':
                print('Head Tail')
        elif whichEnd1 == 'Middle':
            if whichEnd2 == 'Head':
                print('Middle Head')
            elif whichEnd2 == 'Middle':
                print('Middle Middle')
                continue
            elif whichEnd2 == 'Tail':
                print('Middle Tail')
        elif whichEnd1 == 'Tail':
            if whichEnd2 == 'Head':
                print('Tail Head')
            elif whichEnd2 == 'Middle':
                print('Tail Middle')
            elif whichEnd2 == 'Tail':
                print('Tail Tail')


def findUnaligned(hybridfa, seqData, aligndf, paths_c, novatigs, canutigs):
    #INVESTIGATING UNALIGNED CONTIGS FROM SUPERNOVA AND CANU
    param = Parameters()
    hybridData = io.read_fasta(param.hybridFasta)
    seqData.update(hybridData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    hybLenData = {x : len(hybridData[x]) for x in hybridData.keys()}
    
    aligncanu = []
    alignnova = []
    for row in aligndf.iterrows():
        canutig = row[1]['rid']
        novatig = row[1]['qid']
        if str(canutig) not in aligncanu:
            aligncanu.append(str(canutig))
        if str(novatig) not in alignnova:
            alignnova.append(str(novatig))
            
    pathcanu = []
    pathnova = []
    for i in range(0, len(paths_c)):
        for j in range(0, len(paths_c[i])):
            canutig = paths_c[i][j].rid
            novatig = paths_c[i][j].qid
            if str(canutig) not in pathcanu:
                pathcanu.append(str(canutig))
            if str(novatig) not in pathnova:
                pathnova.append(str(novatig))
            
    lonelyNovaAlign = []
    lonelyCanuAlign = []
    lonelyNovaPath = []
    lonelyCanuPath = []
    for tig in novatigs:
        if tig not in alignnova:
            lonelyNovaAlign.append(tig)
        if tig not in pathnova:
            lonelyNovaPath.append(tig)
            
    for tig in canutigs:
        if tig not in aligncanu:
            lonelyCanuAlign.append(tig)
        if tig not in pathcanu:
            lonelyCanuPath.append(tig)
    
    lonelyCanuPathOnly = list(set(lonelyCanuPath).difference(set(lonelyCanuAlign)))
    lonelyNovaPathOnly = list(set(lonelyNovaPath).difference(set(lonelyNovaAlign)))
    
    lens_nova_align = pd.DataFrame(columns=['tig', 'len'])
    lens_nova_path = pd.DataFrame(columns=['tig', 'len'])
    lens_canu_align = pd.DataFrame(columns=['tig', 'len'])
    lens_canu_path = pd.DataFrame(columns=['tig', 'len'])
    ind = 0
    for tig in lonelyNovaAlign:
        lens_nova_align.loc[ind] = tig, lengthData[tig]
        ind += 1
    ind = 0
    for tig in lonelyNovaPathOnly:
        lens_nova_path.loc[ind] = tig, lengthData[tig]
        ind += 1
    ind = 0
    for tig in lonelyCanuAlign:
        lens_canu_align.loc[ind] = tig, lengthData[tig]
        ind += 1
    ind = 0
    for tig in lonelyCanuPathOnly:
        lens_canu_path.loc[ind] = tig, lengthData[tig]
        ind += 1
        
    hybtigList = [26, 31, 79, 95, 98, 99, 108, 134, 244, 256, 290, 327, 382]
    hybtigDict = {}
    for i in hybtigList:
        currentTig = 'hybrid_'+str(i)
        hybtigDict[str(i)+'_region1'] = region_extract(seqData, currentTig, 1, 50000)
        if lengthData[currentTig] > 60000:
            if lengthData[currentTig] < 150000:
                hybtigDict[str(i)+'_region2'] = region_extract(seqData, currentTig, lengthData[currentTig]-50000, lengthData[currentTig])
            else:
                hybtigDict[str(i)+'_region2'] = region_extract(seqData, currentTig, lengthData[currentTig]/2, (lengthData[currentTig]/2)+50000)
                hybtigDict[str(i)+'_region3'] = region_extract(seqData, currentTig, lengthData[currentTig]-50000, lengthData[currentTig])
    
    
    #search all contigs for a sequence
    targetTigs = []
    targetTigs2 = []
    for tig in seqData:
        #print(tig)
        if 'CCAGAGGATTCTTTGCTGTGGGAGGCTGCCCTAGCAATGCTAGGTGTTTCGTTTGACCTCTAAATTT' in seqData[tig]:
            targetTigs.append(tig)
        if 'AAATTTAGAGGTCAAACGAAACACCTAGCATTGCTAGGGCAGCCTCCCACAGCAAAGAATCCTCTGG' in seqData[tig]: #reverse comp
            targetTigs2.append(tig)
