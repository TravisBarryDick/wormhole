#!/usr/bin/bash


dataset=$1
n=$2
numPart=`echo 2^$2 | bc`


#Split train


# Configuration stuff

fspec="./learn/data/$dataset/train/0"
fspect="./learn/data/$dataset/test/0"
outpath="./learn/data/$dataset/random/$numPart"
confIn="./learn/data/$dataset/conf/random"
confOut="./learn/data/$dataset/conf/random.temp$numPart"
# Work out lines per file.


total_lines=$(wc -l <${fspec})
((lines_per_file = (total_lines + numPart - 1) / numPart))

# Split the actual file, maintaining lines.

split -d -a 4 --lines=${lines_per_file} ${fspec} $outpath/train

# Debug information

echo "Total lines     = ${total_lines}"
echo "Lines  per file = ${lines_per_file}"  


# Work out lines per file.

total_lines_t=$(wc -l <${fspect})
((lines_per_file_t = (total_lines_t + numPart - 1) / numPart))

# Split the actual file, maintaining lines.

split -d -a 4 --lines=${lines_per_file_t} ${fspect} $outpath/test

# Debug information

echo "Total lines     = ${total_lines_t}"
echo "Lines  per file = ${lines_per_file_t}" 

#create conf
for i in `seq -f %04g 0 $((numPart-1))`
do
  echo "Starting $i...."
  #create conf
  echo "train_data = \"$outpath/train$i\" " > $confOut
  echo "val_data = \"$outpath/test$i\" " >> $confOut
  outfilename=$outpath/out$i
  #learn
   tracker/dmlc_mpi.py --hostfile=hostfile/bros$n -n 2 -s 1 bin/linear.dmlc $confOut > $outfilename
  rm $confOut $outpath/train$i $outpath/test$i
done

