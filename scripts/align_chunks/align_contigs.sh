REF=$1
QUERY=$2
OUT=$3

: ${OUT:=./out}
mkdir -p $OUT

#parameters
#---------------------------
split_size=1000
threads=8
step=1500  #20000
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

for i in `seq ${start} ${step} ${max}`;
do

    to=$((i + step - 1))
    to=$(($max<$to?$max:$to))
    echo $i to $to
 
    echo dir=${OUT},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads}
    if [ ! -e ${OUT}/alignments/blastn.vs_ref${i}.${to}.out ]
    then
       #echo "file not found."
       qsub -v dir=${OUT},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads} \
	-l nodes=1:ppn=$threads -l mem=${blast_mem}g -l vmem=${blast_mem}g \
        -l walltime=$walltime -o $jobout -e $jobout -d `pwd` \
        -N blast_${i}_${to} sh_runblast
    else
	echo "file exists:${i}-${to}"

    fi  	 

done


