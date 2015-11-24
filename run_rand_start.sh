dataset=$1
n=$2
start=$3
numPart=`echo 2^$2 | bc`


#Split train


# Configuration stuff

fspec="./learn/data/$dataset/train/0"
fspect="./learn/data/$dataset/test/0"
outpath="./learn/data/$dataset/random/$numPart"
confIn="./learn/data/$dataset/conf/random"
confOut="./learn/data/$dataset/conf/random.temp$numPart"


for i in `seq -f %04g $start $((numPart-1))`
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

