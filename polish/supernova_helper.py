import pandas as pd

ID="id"
LEFT="left"
RIGHT="right"

def create_graph_df(fastaIds): 
    
    def parse_10x_name(fastaId, name):
        info10x = dict()
        info10x[ID] = str(fastaId)
        info10x[LEFT] = int(name.split("left=", 1)[1].split(" ", 1)[0])
        info10x[RIGHT] = int(name.split("right=", 1)[1].split(" ", 1)[0])
        return (info10x)
    
    info10x = [parse_10x_name(fid, fastaIds[fid]) for fid in fastaIds]
    graph = pd.DataFrame(info10x)
    return graph

def find_left_contig(tigId, graph):
    leftVertex = int(graph[graph[ID] == tigId][LEFT])
    leftTigs = graph[graph[RIGHT] == leftVertex]
    return list(leftTigs[ID])
    
def find_right_contig(tigId, graph):
    rightVertex = int(graph[graph[ID] == tigId][RIGHT])
    rightTigs = graph[graph[LEFT] == rightVertex]
    return list(rightTigs[ID])
    
def find_haplotig(tigId, graph):
    leftVertex = int(graph[graph[ID] == tigId][LEFT])
    rightVertex = int(graph[graph[ID] == tigId][RIGHT])

    haplotigs =  graph[(graph[LEFT] == leftVertex) & \
                       (graph[RIGHT] == rightVertex)]
    haplotigs = haplotigs[haplotigs[ID] != tigId]
    
    ids = list(haplotigs[ID])
    if len(ids) > 1: print("Warning: more than 1 Supernova haplotig found!")
    if len(ids) < 1: return None

    return ids[0]

