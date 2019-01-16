canu=$1
supernova=$2

out=$3
: ${out:=./out}

mkdir -p $out

#parameters
#---------------------------
split_size=1000
threads=8
step=500  #20000
blast_mem=32  #gb
walltime=7:59:00
jobout=$out/jobout
mkdir -p $jobout

db_dir=$out/blast-db
db=canu_db


#make blast database of canu
#---------------------------
module load blast+/2.7.1
mkdir -p $db_dir
#makeblastdb -in $canu -dbtype nucl -out $db_dir/$db


#split up supernova
#---------------------------
LEN=$split_size
IN=$supernova
mkdir -p $out/split
split_file=$out/split/supernova.${split_size}bp.fa

#cat $IN | perl -nae 'if($_=~/^>/){ if($.!=1){ print "\n"; } print $_; } else{ chomp; print $_; }' | \
#perl -nae 'BEGIN{ $split='$LEN'; } chomp; if($_=~/^>/){ @id=@F; $f=$_; $i=1; $n=0; } else{ if($n==0){ $n=length($_)/1000; @num=split(/\./, $n); if($num[1] ne "00" || $num ne ""){ $n=$num[0]+1; } } for($l=0; $l<=length($_); $l+=$split){ print ">predicted:".substr($id[0],1).":$n:part$i\n".substr($_, $l, $split)."\n"; $i++; } }' > $split_file


#BLAST in parallel
#---------------------------
start=120001 #1
max=`cat $split_file | wc -l`
max=$((max/2))

max=140000

echo $max

for i in `seq ${start} ${step} ${max}`;
do

    to=$((i + step - 1))
    to=$(($max<$to?$max:$to))
    echo $i to $to
 
#    echo dir=${out},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads}
  	 
    qsub -q tcagdenovo -v dir=${out},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads} \
-l nodes=1:ppn=$threads -l mem=${blast_mem}g -l vmem=${blast_mem}g \
-l walltime=$walltime -o $jobout -e $jobout -d `pwd` \
-N blast_${i}_${to} sh_runblast

done
