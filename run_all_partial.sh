dataset=$1

ssh compute-0-17 "cd wormhole/wormhole; nohup time bash run_dddml_partial.sh $dataset 1 > output/$dataset/out/out1 2> output/$dataset/out/2out1" &
#sleep 10m
ssh compute-0-18 "cd wormhole/wormhole; nohup time bash run_dddml_partial.sh $dataset 2 > output/$dataset/out/out2 2> output/$dataset/out/2out2" &
#sleep 20m
ssh compute-0-14 "cd wormhole/wormhole; nohup time bash run_dddml_partial.sh $dataset 3 > output/$dataset/out/out3 2> output/$dataset/out/2out3" &
ssh compute-0-11 "cd wormhole/wormhole; nohup time bash run_dddml_partial.sh $dataset 4 > output/$dataset/out/out4 2> output/$dataset/out/2out4" &
ssh compute-0-12 "cd wormhole/wormhole; nohup time bash run_dddml_partial.sh $dataset 5 > output/$dataset/out/out5 2> output/$dataset/out/2out5" &
ssh compute-0-13 "cd wormhole/wormhole; nohup time bash run_dddml_partial.sh $dataset 6 > output/$dataset/out/out6 2> output/$dataset/out/2out6" &
