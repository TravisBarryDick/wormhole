
#params
conf="./learn/dddml/dddml49.conf"
hostfile="bros_hostfile"
numWorkersDispatch=5
numWorkersDispatchTest=$numWorkersDispatch


#0: clear working directories, make sure folders exist

#1: analyze
learn/dddml/build/analyze_dataset $conf 
#2: sample
learn/dddml/build/datasplit $conf 
#3: cluster
learn/dddml/build/kmeans $conf 
numPart=$?
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
 
