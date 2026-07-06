#!/bin/bash
#SBATCH --job-name=kairos_lr
#SBATCH --output=logs/job_%j_%x.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task 64
#SBATCH --mem=650G
#SBATCH --time=120:00:00
#
#
#
mkdir result/AAAS11	
cp /projects/MA/HGT/wastewater_myloasm/AAAS11/results/AAAS11_assembly.fasta .
snakemake --cores 64 --configfile test-config.yaml --config fasta="AAAS11_assembly.fasta" threads=64 
#mv result/AAAS11* result/AAAS11*
mkdir result/AAI9	
#rm AAAS11_assembly.fasta

cp /projects/MA/HGT/wastewater_myloasm/AAI9/results/AAI9_assembly.fasta .
snakemake --cores 64 --configfile test-config.yaml --config fasta="AAI9_assembly.fasta" threads=64
mv result/AAI9* result/AAI9*

#mkdir AVAS3	
#cp /projects/MA/HGT/wastewater_myloasm/AVAS3/results/AVAS3_assembly.fasta .
#mkdir AVI2	
#cp /projects/MA/HGT/wastewater_myloasm/AVI2/results/AVI2_assembly.fasta .
#mkdir AVI5	
#cp /projects/MA/HGT/wastewater_myloasm/AVI5/results/AVI5_assembly.fasta .
#mkdir EMA16	
#cp /projects/MA/HGT/wastewater_myloasm/EMA16/results/EMA16_assembly.fasta .
#mkdir EMI17	
#cp /projects/MA/HGT/wastewater_myloasm/EMI17/results/EMI17_assembly.fasta .
#mkdir EMI18	
#cp /projects/MA/HGT/wastewater_myloasm/EMI18/results/EMI18_assembly.fasta .
#mkdir FRI15	
#cp /projects/MA/HGT/wastewater_myloasm/FRI15/results/FRI15_assembly.fasta .
#mkdir FRI19	
#cp /projects/MA/HGT/wastewater_myloasm/FRI19/results/FRI19_assembly.fasta .
#mkdir HJ10	
#cp /projects/MA/HGT/wastewater_myloasm/HJ10/results/HJ10_assembly.fasta .
#mkdir HJAS1	
#cp /projects/MA/HGT/wastewater_myloasm/HJAS1/results/HJAS1_assembly.fasta .
#mkdir HJI12	
#cp /projects/MA/HGT/wastewater_myloasm/HJI12/results/HJI12_assembly.fasta .
#mkdir RAS4	
#cp /projects/MA/HGT/wastewater_myloasm/RAS4/results/RAS4_assembly.fasta .
#mkdir RI6	
#cp /projects/MA/HGT/wastewater_myloasm/RI6/results/RI6_assembly.fasta .
#mkdir RI7	
#cp /projects/MA/HGT/wastewater_myloasm/RI7/results/RI7_assembly.fasta .
#mkdir SI13	
#cp /projects/MA/HGT/wastewater_myloasm/SI13/results/SI13_assembly.fasta .
#mkdir SI14	
#cp /projects/MA/HGT/wastewater_myloasm/SI14/results/SI14_assembly.fasta .
#mkdir DAS	
#cp /projects/MA/HGT/wastewater_myloasm/damhusaen_as_rp3/results/damhusaen_as_rp3_assembly.fasta .
#mkdir DI1	
#cp /projects/MA/HGT/wastewater_myloasm/damhusaen_in_rp1/results/damhusaen_in_rp1_assembly.fasta .
#mkdir DI2	
#cp /projects/MA/HGT/wastewater_myloasm/damhusaen_in_rp3/results/damhusaen_in_rp3_assembly.fasta .
#mkdir DIM	
#cp /projects/MA/HGT/wastewater_myloasm/damhusaen_in/results/damhusaen_in_assembly.fasta .

#snakemake --cores 64 --configfile test-config.yaml --config fasta="coded.fasta" threads=64
#snakemake --cores 64 --config fasta="coded.fasta" threads=64
