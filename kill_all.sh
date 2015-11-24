for x in `awk '{if (NR > 8) print $3}' /etc/hosts`; do
	ssh $x "pkill -9 -u vpillutl"
done
