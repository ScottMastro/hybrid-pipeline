REF=$1
QUERY=$2
OUT=$3
PREFIX_NAME="${4:-"query_vs_ref"}"
QUEUE_NAME="${5:-""}"
if [ -z "$QUEUE_NAME" ]; then QUEUE="" ; else QUEUE="-q $QUEUE_NAME"; fi

mkdir -p $OUT
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $SCRIPTDIR

#parameters
#---------------------------
split_size=1000
threads=8
step=20000
blast_mem=32  #gb
walltime=7:59:00
jobout=$OUT/jobout
mkdir -p $jobout

db_dir=$OUT/blast-db
db=ref_db


#make blast database of REF
#---------------------------
module load blast+/2.7.1
mkdir -p $db_dir
makeblastdb -in $REF -dbtype nucl -out $db_dir/$db


#split up QUERY
#---------------------------
LEN=$split_size
IN=$QUERY
mkdir -p $OUT/split
split_file=$OUT/split/query.${split_size}bp.fa

cat $IN | perl -nae 'if($_=~/^>/){ if($.!=1){ print "\n"; } print $_; } else{ chomp; print $_; }' | \
perl -nae 'BEGIN{ $split='$LEN'; } chomp; if($_=~/^>/){ @id=@F; $f=$_; $i=1; $n=0; } else{ if($n==0){ $n=length($_)/1000; @num=split(/\./, $n); if($num[1] ne "00" || $num ne ""){ $n=$num[0]+1; } } for($l=0; $l<=length($_); $l+=$split){ print ">predicted:".substr($id[0],1).":$n:part$i\n".substr($_, $l, $split)."\n"; $i++; } }' > $split_file

#BLAST in parallel
#---------------------------
start=1
max=`cat $split_file | wc -l`
max=$((max/2))

echo $max
JID_LIST=""

for i in `seq ${start} ${step} ${max}`;
do

    to=$((i + step - 1))
    to=$(($max<$to?$max:$to))
    echo $i to $to

    echo dir=${OUT},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads}
    JID=$(qsub -v dir=${OUT},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads} \
     $QUEUE -l nodes=1:ppn=$threads -l mem=${blast_mem}g -l vmem=${blast_mem}g \
     -l walltime=$walltime -o $jobout -e $jobout -d `pwd` \
     -N blast_${i}_${to} ${SCRIPTDIR}/sh_runblast)
     JID_LIST=${JID_LIST}:${JID}


done

JOB="bash ${SCRIPTDIR}/summarize_hits.sh $OUT $PREFIX_NAME"
SUMMARIZE_JID=$(echo $JOB | qsub $QUEUE -W depend=afterok${JID_LIST} -l nodes=1:ppn=1 -l mem=16g -l vmem=16g -l walltime=6:00:00 -o $jobout/ -e $jobout/ -d `pwd` -N summarize_hits "-" | tail -1)


status=`qstat -r $SUMMARIZE_JID | grep $SUMMARIZE_JID | grep [HQR]`
while [ -n "$status" ] # while $status is not empty
   do
      sleep 30 #seconds
      status=`qstat -r $SUMMARIZE_JID | grep $SUMMARIZE_JID | grep [HQR]`
   done


#remove intermediate directories
rm -r $DIR/alignments
rm -r $DIR/split
rm -r $DIR/blast-db
