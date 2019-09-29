canu=$1
supernova=$2

out=$3
: ${out:=./out}

mkdir -p $out

#parameters
#---------------------------
split_size=1000
threads=8
step=1500  #20000
blast_mem=32  #gb
walltime=7:59:00
jobout=$out/jobout
mkdir -p $jobout

db_dir=$out/blast-db
db=canu_db

split_file=$out/split/supernova.${split_size}bp.fa


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
 
    echo dir=${out},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads}
    if [ ! -e ${out}/alignments/blastn.vs_canu${i}.${to}.out ]
    then
       #echo "file not found."
       qsub -v dir=${out},db=${db_dir}/$db,in=${split_file},from=${i},to=${to},threads=${threads} \
	-l nodes=1:ppn=$threads -l mem=${blast_mem}g -l vmem=${blast_mem}g \
        -l walltime=$walltime -o $jobout -e $jobout -d `pwd` \
        -N blast_${i}_${to} sh_runblast
    else
	echo "file exists:${i}-${to}"

    fi  	 

done
