
#params
conf="./learn/dddml/dddml.conf"
hostfile="bros_hostfile"
numWorkersDispatch=5


#0: clear working directories, make sure folders exist

#1: analyze
learn/dddml/build/analyze_dataset $conf 
#2: sample
learn/dddml/build/datasplit $conf 
#3: cluster
learn/dddml/build/kmeans $conf 
#4: dispatch
./tracker/dmlc_mpi.py -n $numWorkersDispatch -s 1 --hostfile $hostfile learn/dddml/build/dispatcher $conf

#5 learning:
numPart=9

for i in `seq 0 $numPart`
do
	echo "Learning $i..."
	./learn/dddml/build/learn $conf train $i
done

for i in `seq 0 $numPart`
do
	echo "Predicting $i..."
	echo "./learn/dddml/build/learn $conf test $i"
done


 
