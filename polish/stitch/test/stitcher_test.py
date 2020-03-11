import sys
sys.path.append("..")
sys.path.append("../..")

import os
import copy
import file_handler as io
import parameters
import stitch.build_intervals as builder

def make_contig(tigId):
    p = parameters.get_parameters()

    filename = "test_data_" + tigId + ".txt"
    alndf = io.parse_alignments(filename)

    #save Intervals for plotting
    cList, bList, bTrashList = [], [], []

    mblockList=[]
    length = alndf.iloc[0][1]*1000
    
    for idx, row in alndf.iterrows():
        #chunks
        chunks = builder.construct_chunks(row, p)
        cList.extend(copy.deepcopy(chunks))
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, p)
        bList.extend(copy.deepcopy(blocks))
        bTrashList.extend(copy.deepcopy(trash))

        #megablock
        mblocks = builder.construct_megablocks(blocks, length, p)
        mblockList.extend(mblocks)

    mList = copy.deepcopy(mblockList)
    
    #contig
    contig = builder.construct_contig(mblockList, tigId, length, p)
    
    return [cList, bList, bTrashList, mList, contig]


def create_test_pickles():
    for filename in os.listdir("."):
        if filename.endswith(".txt") and "test_data" in filename:
            tigId = filename.replace(".txt", "").replace("test_data_", "")
            print(tigId)
            results = make_contig(tigId)       
            io.pickle(results, "./" + tigId + ".pickle", overwrite=True)

def load_test_pickle(tigId):
    return io.unpickle("./" + tigId + ".pickle")


def compare_lists(intervals, intervals_):

    for i, interval in enumerate(intervals):
        if not interval == intervals_[i]:
            print("Lists are not equivalent!")
            print("element #" + str(i) + ":")
            #print(intervals_[i])
            print(interval.__repr__() + '\n' + intervals_[i].__repr__())
            #print(interval)
            return False

    return True

def compare(tigId):
    print("Checking " + tigId + "...")

    
    chunks, blocks, trash, megablocks, contig = make_contig(tigId)
    chunks_, blocks_, trash_, megablocks_, contig_ = load_test_pickle(tigId)
    
    print("Checking chunks...")
    ok1 = compare_lists(chunks, chunks_)
    #if not ok: return False

    print("Checking blocks...")
    ok2 = compare_lists(blocks, blocks_)
    #if not ok: return False
    
    #print("Checking blocks (trash)...")
    #ok = compare_lists(trash, trash_)
    #if not ok: return False
    '''
    for i in range(len(trash)): 
        if not trash[i] == trash_[i]:
            print(i)
    '''
    
    print("Checking megablocks...")
    ok3 = compare_lists(megablocks, megablocks_)
    #if not ok: return False
    
    if ok1 and ok2 and ok3:
        print("PASSED")
        return True
        
compare("516")
compare("429")
compare("111828")
compare("379")

'''
compare("131")
compare("90")
compare("390")
compare("256")
compare("223")
'''