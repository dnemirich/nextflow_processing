#!bin/bash

function usage(){
	echo "Usage:$0 [-w PATH TO WORKING DIRECTORY] [-r PATH TO REFERENCE] [-o PATH TO OUTPUT DIRECTORY]"
}


while getopts "w:r:o:h" opt; do
	case "$opt" in
	w) wd="$OPTARG";;
	r) ref="$OPTARG";;
	o) outdir="$OPTARG";;
	h) usage && exit 0;;
	esac
done


TAGS=$(ls $wd | grep 'bam' | cut -d'.' -f1,1 |sort |uniq)

for TAG in $TAGS;do

picard MergeBamAlignment \
      ALIGNED=$wd/$TAG.sorted.bam  \
      UNMAPPED=$wd/$TAG.unaligned.sorted.bam \
      O=$outdir/$TAG.merged.bam \
      R=$ref \
      SO=coordinate
done

