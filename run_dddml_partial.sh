
#params
dataset=$1
expt=$2
conf="./learn/data/$dataset/conf/dddml$expt.conf"
hostfile="hostfile/bros$expt"
numWorkersDispatch=1
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
./tracker/dmlc_local.py -n $numWorkersDispatch -s 1 learn/dddml/build/dispatcher $conf train
./tracker/dmlc_local.py -n $numWorkersDispatchTest -s 1 learn/dddml/build/dispatcher $conf test


learn/dddml/build/trial $dataset $expt $numPart train
learn/dddml/build/trial $dataset $expt $numPart test


cd local_global/julia-local-global/
julia run_sgd.jl $dataset $expt $numPart 20
cd ../..
exit 0

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
