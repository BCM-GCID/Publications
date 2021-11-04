#!/bin/bash
#PBS -N Busco0
#PBS -o Busco0.out
#PBS -e Busco0.err
#PBS -l nodes=1:ppn=10,mem=20gb
#PBS -q analysis
#PBS -A proj**
source activate py36

busco -f -l eukaryota_odb10 -c 10 -i crypto_BCM2021_v2_min100k.fasta -o crypto_BCM2021_v2_min100k.fasta.busco --config config1 --mode geno
