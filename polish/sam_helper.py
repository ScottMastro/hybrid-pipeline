
def split_haplotype(alignments):
    
    hp1, hp2, unphased = dict(), dict(), []
    
    
    for x in alignments:
        if not x.has_tag("PS") or not x.has_tag("HP"):
            unphased.append(x)
            continue
        
        PS, HP = x.get_tag("PS"), x.get_tag("HP")
            
        if PS not in hp1:
            hp1[PS] = []
            hp2[PS] = []
            
        if HP == 1:
            hp1[PS].append(x)
        elif HP == 2:
            hp2[PS].append(x)
        else:
            print("Unknown HP tag value (" + str(HP) + ") for read " + print(x.to_dict()["name"]))

            
    return (hp1, hp2, unphased)



