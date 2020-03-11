from glob import glob
import pandas as pd
import paf_helper as paf
import gfapy
import re

extension="_pilon_pilon"

def parse_unitigs(canuDir):
    files = glob(canuDir + "/*unitigs.bed")
    if len(files) != 1:
        print("Issue when looking for unitigs.bed file in " + canuDir)
        return None
    
    unitigFile = files[0]
    return pd.read_csv(unitigFile, sep='\t', header=None, names=["chrom", "chromStart", "chromEnd", "name", "score", "strand"])

def parse_gfa(canuDir, contig=False):
    #note: contig gfa does not exist in newer versions of canu

    n = "contigs" if contig else "unitigs" 
    files = glob(canuDir + "/*" + n + ".gfa")
    if len(files) != 1:
        print("Issue when looking for unitigs.gfa file in " + canuDir)
        return None
    
    graphFile = files[0]
    gfa = gfapy.Gfa()

    #canu is dumb and sometimes outputs redundant links????
    with open(graphFile) as graph:
        line = graph.readline()
        uniqueLines = dict()
        while line:
           
            if line[0] == 'L':
                lineSplit = line.split('\t')[1:4]
                if lineSplit[1] == "-":
                    lineSplit[0], lineSplit[2] =  lineSplit[2], lineSplit[0]
                    lineSplit[1] = "+"
                lineKey = "".join(lineSplit)
                if lineKey in uniqueLines: 
                    line = graph.readline()
                    continue
                uniqueLines[lineKey] = True
            
            gfa.add_line(line)
            line = graph.readline()

    graph.close()

    return gfa

def parse_read_map(canuDir, contig=False):

    n = "contigs" if contig else "unitigs" 
    files = glob(canuDir + "/*" + n + ".layout.readToTig")
    if len(files) != 1:
        print("Issue when looking for unitigs.gfa file in " + canuDir)
        return None
    
    readMapFile = files[0]
    
    return pd.read_csv(readMapFile, sep='\t')

def ctg_name(tigId):
    return "ctg" + re.sub("\D", "", str(tigId)).zfill(8).replace(extension, "")
def utg_name(tigId):
    return "utg" + re.sub("\D", "", str(tigId)).zfill(8).replace(extension, "")
def tig_name(tigId, ext=True):
    return "tig" + re.sub("\D", "", str(tigId)).zfill(8) + \
        (extension if ext else "")
def tig_num(tigId):
    return int(re.sub("\D", "", str(tigId).replace(extension, "")))


def map_unitig_contig_custom(canuDir):

    uMap = parse_read_map(canuDir, contig=False)
    cMap = parse_read_map(canuDir, contig=True)

    readMap = cMap.merge(uMap, on="#readID", how='right', suffixes=("_ctg", "_utg"))    
    readMap["minpos_ctg"] = readMap[["bgn_ctg","end_ctg"]].min(axis=1)
    readMap["maxpos_ctg"] = readMap[["bgn_ctg","end_ctg"]].max(axis=1)

    tigMapMin = readMap.groupby(["tigID_utg", "tigID_ctg"])["minpos_ctg"].min()
    tigMap = tigMapMin.index.to_frame()
    tigMap = tigMap.reset_index(drop=True)
    tigMap.columns = ["name", "chrom"]
    tigMap["start"] = tigMapMin.reset_index(drop=True).astype("int32")
    tigMapMax = readMap.groupby(["tigID_utg", "tigID_ctg"])["maxpos_ctg"].max()
    tigMapMax.columns = ["name", "chrom", "end"]
    tigMap["end"] = tigMapMax.reset_index(drop=True).astype("int32")
    tigMap = tigMap[["chrom", "start", "end", "name"]]
    
    tigMap["chrom"] = tigMap["chrom"].apply(ctg_name)
    tigMap["name"] = tigMap["name"].apply(utg_name)
    return tigMap

# bug in gfapy code so I fixed it here...
def is_cut_segment(gfa, segment):
    """Does the removal of a segment split a connected component?
    Note:
      only dovetail overlaps are considered as connections
    Parameters:
      segment (str, Line) : a segment name or instance
    Returns:
       bool
    """
    if isinstance(segment, str):
      segment = gfa.try_get_segment(segment)
    if segment._connectivity() in [(0,0),(0,1),(1,0)]:
      return False
    start_points = list()
    for et in ["L", "R"]:
      for l in segment.dovetails_of_end(et):
        end = l.other_end(gfapy.SegmentEnd(segment.name, et)).inverted()
        if end not in start_points:
            start_points.append(end)
        
    cc = []
    for start_point in start_points:
      cc.append(set())
      visited = set()
      visited.add(segment.name)
      __traverse_component(start_point, cc[-1], visited)
    return any(c != cc[0] for c in cc)

def __traverse_component(segment_end, c, visited):
    s = segment_end.segment
    assert(isinstance(s, gfapy.Line))
    for l in s.dovetails_of_end(segment_end.end_type):
      oe = l.other_end(segment_end)
      sn = oe.name
      s = oe.segment
      if sn in visited:
        continue
      visited.add(sn)
      c.add(s)
      for e in ["L","R"]:
        __traverse_component(gfapy.SegmentEnd(s, e), c, visited)




class CanuGraph:
    def __init__(self, canuDir): 
        self.unitigBed = parse_unitigs(canuDir)
        self.unitigBedSecondary = map_unitig_contig_custom(canuDir)
        self.gfa = parse_gfa(canuDir, contig=False)

    def get_all_unitigs(self):
        return list(set(self.utigAlign["qid"]))

    def get_unitig_neighbours(self, unitigId):
        
        leftTigs = list(set([utg_name(x.name) for x in self.gfa.segment(tig_name(unitigId, ext=False)).neighbours_L]))
        rightTigs = list(set([utg_name(x.name) for x in self.gfa.segment(tig_name(unitigId, ext=False)).neighbours_R]))
        return (leftTigs, rightTigs)
    
    def get_connected_contigs(self, contigId):
        ctg = ctg_name(contigId)
        bed=self.unitigBedSecondary

        allUnitigs = bed[bed["chrom"] == ctg]
        allUnitigIds = set(allUnitigs["name"])
        
        neighbours = set()
        for uid in allUnitigIds:
            l,r = self.get_unitig_neighbours(uid)
            newNeighbours = set([x for x in (l+r) if x not in allUnitigIds])
            neighbours |= newNeighbours

        containedNeighbours = set()
        complexNeighbours = set()
        for neighbour in neighbours:
            l,r =self. get_unitig_neighbours(neighbour)
            newNeighbours = set([x for x in (l+r) if x not in allUnitigIds])
            if len(newNeighbours) > 0:
                complexNeighbours.add(neighbour)
            else:
                containedNeighbours.add(neighbour)
        
        containedContigs = set(bed[bed["name"].isin(containedNeighbours)]["chrom"])
        containedContigIds = [tig_name(x) for x in containedContigs]
        
        complexContigs = set(bed[bed["name"].isin(complexNeighbours)]["chrom"])
        complexContigIds = [tig_name(x) for x in complexContigs]
        
        return(containedContigIds, complexContigIds)
                    

    def get_connected_contigs2(self, contigId):
        
        bed=refGraph.unitigBedSecondary
        refgfa=refGraph.gfa
        
        ctg = ctg_name(contigId)
        unitigs = bed[bed["chrom"] == ctg]
        unitigIds = list(unitigs["name"])
        seg = refgfa.segment(tig_name(unitigIds[0], ext=False))
        connectedSegments = refgfa.segment_connected_component(seg)
        connectedSegmentsIds = {utg_name(seg.name) for seg in connectedSegments}

        contigs = bed[bed["name"].isin(connectedSegmentsIds)]
        contigMap = {x["name"] : x["chrom"] for _,x in contigs.iterrows()}
        contigIds = set(contigs["chrom"])
        
        contigsStarts = contigs.groupby(["chrom"]).min()       
        contigsEnds = contigs.groupby(["chrom"]).max()       

        starts = { x["name"] : x.name for _,x in contigsStarts.iterrows() }
        ends = { x["name"] : x.name for _,x in contigsEnds.iterrows() }












        leftTigs = { cid : set() for cid in contigIds }
        rightTigs = { cid : set() for cid in contigIds }

        for uid in connectedSegmentsIds: #[seg.name for seg in connectedSegments]:
            l,r = refGraph.get_unitig_neighbours(uid)
            cid = contigMap[uid]
            for uidL in l:
                cidL = contigMap[uidL]
                if cid != cidL: leftTigs[cid].add(cidL)
            for uidR in r:
                cidR = contigMap[uidR]
                if cid != cidR: rightTigs[cid].add(cidR)


        

        contigs = bed[bed["name"].isin(connectedSegmentsIds)]





        for uid in connectedSegmentsIds: #[seg.name for seg in connectedSegments]:
            l,r = refGraph.get_unitig_neighbours(uid)
            cid = contigMap[uid]
            for uidL in l:
                cidL = contigMap[uidL]
                if cid != cidL: leftTigs[cid].add(cidL)
            for uidR in r:
                cidR = contigMap[uidR]
                if cid != cidR: rightTigs[cid].add(cidR)
   
        containedTigs = dict()
        # figure out easy ones
        for cid in contigIds:
            l, r = leftTigs[cid], rightTigs[cid] 
            if len(l) == 1 and len(r) == 1 and l == r:
                containedTigs[cid] = list(l)[0]
               
        containedTigList = [c for c in containedTigs]
        for contained in containedTigList:
            leftTigs.pop(contained)
            rightTigs.pop(contained)
            container = containedTigs[contained]
            leftTigs[container].remove(contained)
            rightTigs[container].remove(contained)
            
        
        starts = [cid for cid in leftTigs if len(leftTigs[cid]) == 0]
        ends = [cid for cid in rightTigs if len(rightTigs[cid]) == 0]
        middle = [seg for seg in connectedSegments if (seg not in starts and seg not in ends)]


            
            
            
                
                
        
         
        unitigIds1 = list(unitigs1["name"])

        for startSegment in starts:
            uid = utg_name(startSegment.name)
            traversed = { uid } 
            
            while True:
                neighbours = set()
                
                
                
                
                
                newNeighbours = set([x for x in (l+r) if x not in allUnitigIds])
                neighbours |= newNeighbours

            
            
            for uid in allUnitigIds:
    
                containedNeighbours = set()
                complexNeighbours = set()
                for neighbour in neighbours:
                    l,r =self. get_unitig_neighbours(neighbour)
                    newNeighbours = set([x for x in (l+r) if x not in allUnitigIds])
                    if len(newNeighbours) > 0:
                        complexNeighbours.add(neighbour)
                    else:
                        containedNeighbours.add(neighbour)
            
            containedContigs = set(bed[bed["name"].isin(containedNeighbours)]["chrom"])
            containedContigIds = [tig_name(x) for x in containedContigs]
            
            complexContigs = set(bed[bed["name"].isin(complexNeighbours)]["chrom"])
            complexContigIds = [tig_name(x) for x in complexContigs]
            
            return(containedContigIds, complexContigIds)

        
        


        bottleneck = []
        for segment in segments:
            bottleneck.append([segment.name, is_cut_segment(refgfa, segment)])
        
        ctg2 = ctg_name(contigId2)

        
        ctg = ctg_name(contigId)
        bed=self.unitigBedSecondary

        allUnitigs = bed[bed["chrom"] == ctg]
        allUnitigIds = set(allUnitigs["name"])
        
        neighbours = set()
        for uid in allUnitigIds:
            l,r = self.get_unitig_neighbours(uid)
            newNeighbours = set([x for x in (l+r) if x not in allUnitigIds])
            neighbours |= newNeighbours

        containedNeighbours = set()
        complexNeighbours = set()
        for neighbour in neighbours:
            l,r =self. get_unitig_neighbours(neighbour)
            newNeighbours = set([x for x in (l+r) if x not in allUnitigIds])
            if len(newNeighbours) > 0:
                complexNeighbours.add(neighbour)
            else:
                containedNeighbours.add(neighbour)
        
        containedContigs = set(bed[bed["name"].isin(containedNeighbours)]["chrom"])
        containedContigIds = [tig_name(x) for x in containedContigs]
        
        complexContigs = set(bed[bed["name"].isin(complexNeighbours)]["chrom"])
        complexContigIds = [tig_name(x) for x in complexContigs]
        
        return(containedContigIds, complexContigIds)



        
        """
        def build_regions(self):
        
        allUnitigs = self.get_all_unitigs()
        unitigMap = dict(zip(unitigBed["name"], unitigBed["chrom"]))
        unitigLen = dict(zip(unitigBed["name"], (unitigBed["chromEnd"]-unitigBed["chromStart"])))
        unitigStartMap = dict(zip(unitigBed["name"], unitigBed["chromStart"]))
        unitigEndMap = dict(zip(unitigBed["name"], unitigBed["chromEnd"]))

        alignGroups = utigAlign.groupby(["qid"])

        for unitigId in allUnitigs:
            utg = utg_name(unitigId)
            
            #check if bed file is helpful:        
            if utg in unitigMap:
                ctg = unitigMap[utg]
            else:
                print("missing")
                continue
            
            
            
            
            
            
            



            contigId = tig_name(ctg)
        
            group = alignGroups.get_group(unitigId)
            group = group[group["rid"] == contigId]
            if len(group) == 0:
                print("empty")
                continue
            
            start = min(min(group["rstart"]),  min(group["rend"]))
            end = max(max(group["rstart"]),  max(group["rend"]))
            
            ustart = min(min(group["qstart"]),  min(group["qend"]))
            uend = max(max(group["qstart"]),  max(group["qend"]))
            uregion = SimpleRegion(unitigId, ustart, uend)

            region = SimpleRegion(contigId, start, end)
            lengthRatio = len(region)/unitigLen[utg]
            ulengthRatio = len(uregion)/unitigLen[utg]

            if lengthRatio < 0.95 or lengthRatio > 1.05:
                print(lengthRatio, ulengthRatio)
                input()
        
            alignments = self.utigAlign[(self.utigAlign["rid"] == ctg) & (self.utigAlign["qid"] == utg)]
            
            start = min(min(alignments["rstart"]),  min(alignments["rend"]))
            end = max(max(alignments["rstart"]),  max(alignments["rend"]))
    
            SimpleRegion(contigId, start, end)
            
            
            
            
            
            
                
            contigIds = unitigBed[unitigBed["name"] == utg]
            
            [tig_name(x.chrom) for x in unitigBed.filter(lambda b: b.name == utg)]
            
            if len(contigIds) > 1:
                print("wait") ; input()
            
                    
        
        return list(set(self.utigAlign["qid"]))
        
        
        
            def contig_to_unitigs(self, contigId, idsOnly=True):
        ctg = ctg_name(contigId)
        
        if idsOnly:  return [x.name for x in self.unitigMap.filter(lambda b: b.chrom == ctg)]
        return pybedtools.BedTool([x for x in self.unitigMap.filter(lambda b: b.chrom == ctg)])
    
    def unitig_to_contig(self, unitigId, idsOnly=True):
        utg = utg_name(unitigId)
        
        if idsOnly:
            match = [tig_name(x.chrom) for x in self.unitigMap.filter(lambda b: b.name == utg)]
            if len(match) == 0:
                return[tig_name(x.chrom) for x in self.unitigMapSecondary.filter(lambda b: b.name == utg)]
        else:
            match = pybedtools.BedTool([x for x in self.unitigMap.filter(lambda b: b.name == utg)])
            if len(match) == 0:
                match = pybedtools.BedTool([x for x in self.unitigMapSecondary.filter(lambda b: b.name == utg)])            
        
        return match

    def get_closest_unitig(self, contigId, pos):
        
        '''
        unitigsBed = self.contig_to_unitigs(contigId, idsOnly=False)
        
        point = pybedtools.BedTool(" ".join([contigId, str(pos), str(pos)]), from_string=True)
        closest = unitigsBed.closest(point)
        
        if len(closest) < 1: return None
        return closest[0].name
        '''
        tig = tig_name(contigId)
        
        alignments = self.utigAlign[(self.utigAlign["rid"] == tig)]
        uids = set(alignments["qid"])
        
        regions = [self.unitig_alignment_region(contigId, uid) for uid in uids]
        
        alignments["sdist"] = self.utigAlign["rstart"] -pos


    def unitig_alignment_region(self, contigId, unitigId):
        utg = tig_name(unitigId, ext=False)
        ctg = tig_name(contigId)
    
        alignments = self.utigAlign[(self.utigAlign["rid"] == ctg) & (self.utigAlign["qid"] == utg)]
        
        start = min(min(alignments["rstart"]),  min(alignments["rend"]))
        end = max(max(alignments["rstart"]),  max(alignments["rend"]))

        return SimpleRegion(contigId, start, end)
        
    def unitig_pos(self, contigId, unitigId):
        utg = utg_name(unitigId)
        ctg = ctg_name(contigId)
    
        match = [x for x in self.unitigMap.filter(lambda b: b.chrom == ctg and b.name == utg)]
        
        if len(match) < 1:
            match = [x for x in self.unitigMapSecondary.filter(lambda b: b.chrom == ctg and b.name == utg)]
        
        if len(match) < 1: return (None, None)
        if len(match) > 1: 
            print("Warning: multiple matches for contig-unitig pair " + ctg + " " + utg)
            
        return (int(match[0][1]), int(match[0][2]))        
        
    """