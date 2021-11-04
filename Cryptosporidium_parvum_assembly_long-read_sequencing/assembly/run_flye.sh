#!/bin/bash
#PBS -N Flye
#PBS -o flye.out
#PBS -e flye.err
#PBS -l nodes=1:ppn=15,mem=90gb
#PBS -l walltime=44:00:00
#PBS -q analysis
#PBS -A proj**

flye --meta --nano-raw pass_all.fastq -g 9m  -o output_flye_pass_meta -t 15 -i 2
