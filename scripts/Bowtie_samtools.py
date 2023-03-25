#!/usr/bin/env python

import argparse
import subprocess
import os

# initiate argument parsing
parser = argparse.ArgumentParser()

# add flags
parser.add_argument('--sample_id', required=True, help='Sample identifier')
parser.add_argument('--read', required=True, help='Input reads file')

# add flag input to variable
args = parser.parse_args() 

reference_index = 'MG1655'

# bowtie2 alignment
bowtie_command = f'bowtie2 --no-unal -x {reference_index} -U {args.read} -p 10 -S {args.sample_id}.sam 2> {args.sample_id}_.txt'
try:
    subprocess.run(bowtie_command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    print(f'Error: bowtie2 command failed with error code {e.returncode}')
    exit(e.returncode)

# samtools view
samtools_view_command = f'samtools view --threads 10 {args.sample_id}.sam -o {args.sample_id}.bam'
try:
    subprocess.run(samtools_view_command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    print(f'Error: samtools view command failed with error code {e.returncode}')
    exit(e.returncode)

# samtools sort
samtools_sort_command = f'samtools sort {args.sample_id}.bam -o {args.sample_id}.sorted.bam'
try:
    subprocess.run(samtools_sort_command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    print(f'Error: samtools sort command failed with error code {e.returncode}')
    exit(e.returncode)

# remove intermediate files
os.remove(f'{args.sample_id}.bam')
os.remove(f'{args.sample_id}.sam')