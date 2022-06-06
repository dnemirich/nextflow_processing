#!bin/bash

function usage(){
	echo "Usage:$0 [-f PATH TO FORWARD READS] [-r PATH TO REVERSE READS] [-T THRESHOLD FOR FASTQ-FILTERER] [-t TILES TO REMOVE] [-a FIRST ADAPTER] [-A SECOND ADAPTER] [-c NUMBER OF CORES] [-o PATH TO OUTPUT DIRECTORY]"
}


while getopts "f:r:T:t:a:A:c:o:h" opt; do
	case "$opt" in
	f) forward="$OPTARG";;
	r) reverse="$OPTARG";;
	T) threshold="$OPTARG";;
	t) tiles="$OPTARG";;
	a) adapter_1="$OPTARG";;
	A) adapter_2="$OPTARG";;
	c) cores="$OPTARG";;
	o) outdir="$OPTARG";;
	h) usage && exit 0;;
	esac
done


forward_out=$(ls $forward | xargs -n 1 basename | sed 's/.fastq.gz//')
reverse_out=$(ls $reverse | xargs -n 1 basename | sed 's/.fastq.gz//')

./Fastq-Filterer/fastq_filterer --i1 $forward --i2 $reverse --threshold $threshold --remove_tiles $tiles --o1 $outdir/${forward_out}_filtered.fastq --o2 $outdir/${reverse_out}_filtered.fastq

cutadapt --discard-trimmed -g $adapter_1 -G $adapter_2 -j $cores -o $outdir/${forward_out}_filtered.noadapt.fastq.gz -p $outdir/${reverse_out}_filtered.noadapt.fastq.gz $outdir/${forward_out}_filtered.fastq $outdir/${reverse_out}_filtered.fastq



