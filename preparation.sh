#!bin/bash

conda create -n melanoma python=3.6 --file packages.txt

conda activate melanoma
git clone https://github.com/EdinburghGenomics/Fastq-Filterer.git
cd Fastq-Filterer
make
cd ..

mkdir -p input
mkdir -p reference
mkdir -p data/funcotator

cd data/funcotator
gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
tar -xf *.tar.gz
cd ..
cd ..

