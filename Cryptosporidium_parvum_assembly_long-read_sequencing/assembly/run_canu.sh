#!/bin/bash
#PBS -N Canu
#PBS -o canu.out
#PBS -e canu.err
#PBS -l nodes=1:ppn=15,mem=90gb
#PBS -l walltime=44:00:00
#PBS -q analysis
#PBS -A proj**


canu -p crypto_pass_ont -d output_pass_canu genomeSize=9m  -nanopore-raw pass_all.fastq useGrid=false
