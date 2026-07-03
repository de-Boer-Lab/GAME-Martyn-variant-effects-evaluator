import json
import pandas as pd
from collections import Counter
import pyfaidx
import kipoiseq
from kipoiseq import Interval
import pandas as pd
import numpy as np
import random
from tqdm import tqdm
import math
import pysam
import glob
import os

def merge_variants_from_tsv(path_to_dataFolder, cell_type):
    fasta_path = "/arc/project/st-cdeboer-1/iluthra/hg38.fa"

    genome = pysam.FastaFile(fasta_path)

    # Find all .tsv files in the folder
    tsv_files = glob.glob(os.path.join(path_to_dataFolder, "*.tsv"))

    # Make sure there are files to read
    if not tsv_files:
        raise ValueError("No .tsv files found in the specified folder.")

    # Read the first file to get the column names and initialize the main DataFrame
    snp_dataFrame = pd.read_csv(tsv_files[0], sep='\t')
    # Read the rest and append
    for file in tsv_files[1:]:
        df = pd.read_csv(file, sep='\t')
        snp_dataFrame = pd.concat([snp_dataFrame, df], ignore_index=True)

    snp_dataFrame.to_csv(path_to_dataFolder + '/all_' + cell_type + '.tsv', sep='\t')

  
#step 1: create a merged TSV for the variant files for each of the cell types
merge_variants_from_tsv('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Engritz_Evaluator/evaluator_data/Jurkat/', "Jurkat")
merge_variants_from_tsv('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Engritz_Evaluator/evaluator_data/THP1/', "THP1")

#step 2: use the spdiToSequence_IL.R script to pull the Reference and alternate sequence for each variant
#ran this locally on IL computer since it only needs to be done once
#conda activate r_SPDI
#Rscript spdiToSequence_IL.R /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/Jurkat/all_Jurkat.tsv 1000 /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/Jurkat/all_Jurkat_sequences.tsv
#Rscript spdiToSequence_IL.R /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/THP1/all_THP1.tsv 1000 /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/THP1/all_THP1_sequences.tsv

