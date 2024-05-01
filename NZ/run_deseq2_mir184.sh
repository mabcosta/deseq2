
reference="scramble"
for condition in $(cut -d, -f2 sample_info.csv | tail -n+2 | grep -v $reference | sort | uniq); do
	echo "-------- RUNNING: $condition and reference $refence -----------"

	Rscript NZ_script.R -c $condition -r $refence
	mkdir -p ~/MarcosCosta/projects/MIR184/DESeq2/results/"$condition"_"$refence" 
	mv mir184_scramble_* ~/MarcosCosta/projects/MIR184/DESeq2/results/"$condition"_"$refence"
	
	echo "-------- FINISHED: $condition and reference $refence -----------"
	echo "----------------------------------------------------------------"
	echo "----------------------------------------------------------------"
done 

reference="mir184"
for condition in $(cut -d, -f2 sample_info.csv | tail -n+2 | grep -v $reference | sort | uniq); do
	echo "-------- RUNNING: $condition and reference $refence -----------"

	Rscript NZ_script.R -c $condition -r $refence
	mkdir -p ~/MarcosCosta/projects/MIR184/DESeq2/results/"$condition"_"$refence" 
	mv mir184_scramble_* ~/MarcosCosta/projects/MIR184/DESeq2/results/"$condition"_"$refence"
	
	echo "-------- FINISHED: $condition and reference $refence -----------"
	echo "----------------------------------------------------------------"
	echo "----------------------------------------------------------------"

done 