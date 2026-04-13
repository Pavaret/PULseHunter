#!/bin/bash
conda create -n pul_hmm_2026 -y
conda activate pul_hmm_2026 
conda install pandas biopython prodigal hmmer mafft mmseqs2 blast -y
