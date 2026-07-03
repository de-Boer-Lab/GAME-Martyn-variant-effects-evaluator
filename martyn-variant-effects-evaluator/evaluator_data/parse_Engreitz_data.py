import pandas as pd
# import pyfaidx
# import kipoiseq
# from kipoiseq import Interval
import pandas as pd
import pysam
import pysam
from Bio.Seq import Seq

# def merge_variants_from_tsv(path_to_dataFolder, cell_type):
    
#     # Find all .tsv files in the folder
#     tsv_files = glob.glob(os.path.join(path_to_dataFolder, "*.tsv"))

#     # Make sure there are files to read
#     if not tsv_files:
#         raise ValueError("No .tsv files found in the specified folder.")

#     # Read the first file to get the column names and initialize the main DataFrame
#     snp_dataFrame = pd.read_csv(tsv_files[0], sep='\t')
#     # Read the rest and append
#     for file in tsv_files[1:]:
#         df = pd.read_csv(file, sep='\t')
#         snp_dataFrame = pd.concat([snp_dataFrame, df], ignore_index=True)

#     snp_dataFrame.to_csv(path_to_dataFolder + '/all_' + cell_type + '.tsv', sep='\t')

  
# #step 1: create a merged TSV for the variant files for each of the cell types
# merge_variants_from_tsv('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Engritz_Evaluator/evaluator_data/Jurkat/', "Jurkat")
# merge_variants_from_tsv('/arc/project/st-cdeboer-1/iluthra/API_genomic_model_evaluation/Engritz_Evaluator/evaluator_data/THP1/', "THP1")

#step 2: use the spdiToSequence_IL.R script to pull the Reference and alternate sequence for each variant
#ran this locally on IL computer since it only needs to be done once
#conda activate r_SPDI
#Rscript spdiToSequence_IL.R /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/Jurkat/all_Jurkat.tsv 1000 /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/Jurkat/all_Jurkat_sequences.tsv
#Rscript spdiToSequence_IL.R /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/THP1/all_THP1.tsv 1000 /Users/ishika/Desktop/API/Engritz_Evaluator/evaluator_data/THP1/all_THP1_sequences.tsv

###Instead of using the SPDI script we are going to center at PPIF/ILR2A Promoters and pull 1MB up/downstream of the promoters
###Variants will be inserted into the sequences
###Prediction_ranges with will over the promoter

#PPIF promoter target region is chr10: 81,107,026-81,107,246 (hg19) 220bp
#IL2RA promoter target region is chr10:6,104,468-6,104,652 (hg19)  184bp
#From paper
PPIF_chrom = "chr10"
PPIF_start = 81107026
PPIF_strand = "positive"
IL2RA_chrom = "chr10"
IL2RA_start = 6104468
IL2RA_strand = "negative"
fasta_path = "/arc/project/st-cdeboer-1/Genomes/human/GRChg37/hg19.fa"
genome = pysam.FastaFile(fasta_path)

PPIF_reference_seq = genome.fetch(PPIF_chrom, PPIF_start-1000000,  PPIF_start+1000000)
PPIF_reference_seq  = PPIF_reference_seq.upper()
print(len(PPIF_reference_seq))
IL2RA_reference_seq = genome.fetch(IL2RA_chrom, IL2RA_start-1000000,  IL2RA_start+1000000)
###Need to take the reverse complement of the IL2RA reference sequence
IL2RA_reference_seq_negative_strand = Seq(IL2RA_reference_seq).reverse_complement()
IL2RA_reference_seq  = IL2RA_reference_seq.upper()
IL2RA_reference_seq_negative_strand = IL2RA_reference_seq_negative_strand.upper()
print(len(IL2RA_reference_seq))

# Build a dictionary where keys are gene names and values are sequences
seq_dict = {
    "PPIF": PPIF_reference_seq,
    "IL2RA": IL2RA_reference_seq_negative_strand
}

# Convert to DataFrame
df = pd.DataFrame.from_dict(seq_dict, orient="index", columns=["sequence"])
df['gene_length'] = [200, 184]
df['strand'] = ["positive", "negative"]
# Optional: reset index to have 'gene' as a column instead of index
df.index.name = "gene"
df.reset_index(inplace=True)
df = df[["gene", "gene_length", "strand", "sequence"]]
# Write to CSV
df.to_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Engreitz_Evaluator_final/evaluator_data/gene_reference_sequences.csv", index=False)

def parse_variant_id(v):
    chrom, pos, alleles = v.split(":")
    ref, alt = alleles.split(">")
    return chrom, int(pos), ref, alt

def insert_variant(seq, backbone_start, pos, ref, alt, strand):
    idx = pos - backbone_start
    seq = list(seq)

    # Safety check
    if seq[idx:idx+len(ref)] != list(ref):
        print(
            f"REF mismatch at {pos}: "
            f"expected {ref}, found {''.join(seq[idx:idx+len(ref)])}"
        )
    seq[idx:idx+len(ref)] = list(alt)
    seq_pulled = "".join(seq)
    #If the gene is on the negative strand you need to take the reverse complement after you put the mutation in
    if strand == "negative":
        seq_pulled = Seq(seq_pulled).reverse_complement()
    return seq_pulled

all_THP1 = pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Engreitz_Evaluator_final/evaluator_data/THP1/all_THP1.tsv", sep = '\t')
all_Jurkat = pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Engreitz_Evaluator_final/evaluator_data/Jurkat/all_Jurkat.tsv", sep = '\t')

variant_sequences = []
prediction_ranges = []
for i, row in all_THP1.iterrows():

    chrom, pos, ref, alt = parse_variant_id(row["VariantID_h19"])
    gene = row["gene_symbol"]

    if gene == "PPIF":
        backbone = PPIF_reference_seq
        backbone_start = PPIF_start -1000000
        expected_chrom = PPIF_chrom
        strand = PPIF_strand
    if gene == "IL2RA":
        backbone = IL2RA_reference_seq 
        backbone_start = IL2RA_start -1000000
        expected_chrom = IL2RA_chrom
        strand = IL2RA_strand
    if chrom != expected_chrom:
        raise ValueError(f"Chrom mismatch for {gene}: {chrom}")

    seq_with_variant = insert_variant(
        backbone,
        backbone_start,
        pos,
        ref,
        alt,
        strand
    )
    #print(len(seq_with_variant))
    #print(seq_with_variant[(pos-backbone_start-10):(pos-backbone_start+10)])
    variant_sequences.append(seq_with_variant)

all_THP1["variant_sequence"] = variant_sequences
print(all_THP1)

all_THP1.to_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Engreitz_Evaluator_final/evaluator_data/THP1/variant_sequences.csv" , sep = '\t')

#Do the same for Jurkat
variant_sequences = []

for i, row in all_Jurkat.iterrows():

    chrom, pos, ref, alt = parse_variant_id(row["VariantID_h19"])
    gene = row["gene_symbol"]

    if gene == "PPIF":
        backbone = PPIF_reference_seq
        backbone_start = PPIF_start -1000000
        expected_chrom = PPIF_chrom
    if gene == "IL2RA":
        backbone = IL2RA_reference_seq 
        backbone_start = IL2RA_start -1000000
        expected_chrom = IL2RA_chrom

    if chrom != expected_chrom:
        raise ValueError(f"Chrom mismatch for {gene}: {chrom}")

    seq_with_variant = insert_variant(
        backbone,
        backbone_start,
        pos,
        ref,
        alt,
        strand
    )
    #print(len(seq_with_variant))
    #print(seq_with_variant[(pos-backbone_start-10):(pos-backbone_start+10)])
    variant_sequences.append(seq_with_variant)

all_Jurkat["variant_sequence"] = variant_sequences
print(all_Jurkat)

all_Jurkat.to_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Engreitz_Evaluator_final/evaluator_data/Jurkat/variant_sequences.csv" , sep = '\t')
