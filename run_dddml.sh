
#params
dataset=$1
expt=$2
conf="./learn/data/$dataset/conf/dddml$expt.conf"
hostfile="hostfile/bros$expt"
numWorkersDispatch=5
numWorkersDispatchTest=$numWorkersDispatch
echo $conf 
#0: clear working directories, make sure folders exist
find ./learn/data/$dataset/experiment$expt -type f -exec rm {} \;


#1: analyze
learn/dddml/build/analyze_dataset $conf 
#2: sample
learn/dddml/build/datasplit $conf 
#3: cluster
learn/dddml/build/kmeans $conf 
numPart=$?
echo "$numPart Clusters..."
echo $numPart > ./learn/data/$dataset/experiment$expt/ncluster
#return value of clustering should be the final number of clusters
#4: dispatch
./tracker/dmlc_mpi.py -n $numWorkersDispatch -s 1 --hostfile $hostfile learn/dddml/build/dispatcher $conf train
./tracker/dmlc_mpi.py -n $numWorkersDispatchTest -s 1 --hostfile $hostfile learn/dddml/build/dispatcher $conf test

#5 learning:

for i in `seq 0 $((numPart-1))`
do
	echo "Learning $i..."
	./learn/dddml/build/learn $conf train $i
done

for i in `seq 0 $((numPart-1))`
do
	echo "Predicting $i..."
	./learn/dddml/build/learn $conf test $i $numWorkersDispatchTest
done

./learn/dddml/build/evaluate_accuracy $conf $numPart $numWorkersDispatchTest

find ./learn/data/$dataset/experiment$expt -type f -exec rm {} \;

echo "DONE $2"
