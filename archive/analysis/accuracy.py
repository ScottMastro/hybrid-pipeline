import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
before = "../input/ngaps_hg38.txt"

after ="../input/nfixed_hg38.txt"
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore",
        "qlen", "slen", "gaps", "qseq", "sseq", "sframe"]

aln_before = pd.read_csv(before, sep="\t", names=cols)
aln_after = pd.read_csv(after, sep="\t", names=cols)



plt.scatter(aln_before["length"], aln_before["qlen"], alpha=0.01)

aln_before["partial"] = (aln_before["length"]) <= 600
aln_before["full"] = abs(aln_before["length"] - aln_before["qlen"]) >= 100
aln_before["same"] = (aln_before["partial"] == 0) & (aln_before["extra"] == 0)
aln_before["fixed"] = aln_before["qseqid"].isin(aln_after["qseqid"])

#partial
partials_before = aln_before[(aln_before["partial"] == True) & (aln_before["fixed"] == True)]
partials_before = partials_before.sort_values(by=['qseqid'])
partials_after = aln_after[aln_after["qseqid"].isin(partials_before["qseqid"])]
partials_after = partials_after.sort_values(by=['qseqid'])

plt.hist(partials_before["pident"], bins=200, alpha = 0.5)
plt.hist(partials_after["pident"], bins=200,  alpha =0.5)

partials_after["diff"] = np.array(partials_after["length"]) - np.array(partials_before["length"])
partials_after["success"] = partials_after["length"] > 1000
partials_after["unclear"] = (partials_after["length"] < 1000) & (partials_after["diff"] > 100)
partials_after["fail"] = (partials_after["diff"] <= 100)

print("SUCCESS: " + str(partials_after["success"].sum()) + 
      " (" + str(round(100.0*partials_after["success"].sum()/len(partials_after) ,4)) + "%)" )
print("UNCLEAR: " + str(partials_after["unclear"].sum()) + 
      " (" + str(round(100.0*partials_after["unclear"].sum()/len(partials_after) ,4)) + "%)" )
print("FAIL: " + str(partials_after["fail"].sum()) + 
      " (" + str(round(100.0*partials_after["fail"].sum()/len(partials_after) ,4)) + "%)" )

##
plt.scatter(partials_after[partials_after["fail"] ]["length"], 
            partials_after[partials_after["fail"] ]["qlen"], alpha=0.2)
plt.xlabel('Alignment Length', fontsize=18)
plt.ylabel('Query Length', fontsize=16)
##

##
plt.hist(partials_after[partials_after["fail"] ]["pident"], bins=30)
plt.xlabel('Alignment % Identity', fontsize=18)
##

partials_after["unclear_ok"] = ((np.array(partials_after["qlen"]) -
              np.array(partials_after["length"])) < 100) & (partials_after["unclear"])

##
plt.scatter(partials_after[~partials_after["unclear_ok"] &partials_after["unclear"] ]["length"], 
            partials_after[~partials_after["unclear_ok"] &partials_after["unclear"]]["qlen"], alpha=0.2)
plt.xlabel('Alignment Length', fontsize=18)
plt.ylabel('Query Length', fontsize=16)
##


#

aln_after["diff"] = aln_after["qlen"] - aln_after["length"]

full_before = aln_before[(aln_before["full"] == True) & (aln_before["fixed"] == True)]
full_before = full_before.sort_values(by=['qseqid'])
full_after = aln_after[aln_after["qseqid"].isin(partials_before["qseqid"])]
full_after = partials_after.sort_values(by=['qseqid'])

plt.hist(aln_after[aln_after["diff"] < 1000]["diff"], bins=50)
plt.xlabel("(query length) - (align length)", fontsize=18)

plt.hist(aln_after["pident"], bins=30)
plt.xlabel('Alignment % Identity', fontsize=18)

plt.scatter(aln_after["length"], aln_after["qlen"], c=aln_after["pident"], alpha=0.5)











def meanSd(aln_before, aln_after, column):
    values_before = []
    values_after = []
     
    for index, row_after in aln_after.iterrows():
        row_id = row_after["qseqid"]
        row_before = aln_before[aln_before["qseqid"]==row_id].iloc[0]
       
        values_before.append(row_before[column])
        values_after.append(row_after[column])
               
    print(column + "=" + str(np.mean(values_before)) + " ~ " + str(np.std(values_before)))
    print(column + "=" + str(np.mean(values_after)) + " ~ " + str(np.std(values_after)))
    plt.hist(values_before)
    plt.hist(values_after)



def count(aln_before, aln_after, column, cmpr="greater"):
    counter= 0
    n=0
     
    for index, row_after in aln_after.iterrows():
        row_id = row_after["qseqid"]
        row_before = aln_before[aln_before["qseqid"]==row_id].iloc[0]
    
        if cmpr == "greater" and row_before[column] < row_after[column] :
            counter=counter+1
        elif cmpr == "less" and row_before[column] > row_after[column] :
            counter=counter+1

        n=n+1
               
    if cmpr == "greater":
        print(column + " increased for " + str(counter) + "/" + str(n) + " (" + str(100*counter/n) + "%)")
    elif cmpr == "less":
        print(column + " increased for " + str(counter) + "/" + str(n) + " (" + str(100*counter/n) + "%)")



length_increase = []
lens_before = []

deltas = []
pid_change = []
shrink_counter= 0
n=0
for index, row_after in aln_after.iterrows():
    n=n+1
    row_id = row_after["qseqid"]
    row_before = aln_before[aln_before["qseqid"]==row_id].iloc[0]

    before_len = row_before["qlen"]
    after_len = row_after["qlen"]
    
    delta = after_len - before_len
    deltas.append(delta)
    if(delta < 0):
       shrink_counter= shrink_counter+1
       pid_change.append(row_after["pident"] - row_before["pident"])
    else:
       length_increase.append(row_after["length"] - delta - row_before["length"])
       lens_before.append(row_before["length"])
           
plt.scatter(length_increase, lens_before, alpha=0.1)
plt.hist(pid_change)
