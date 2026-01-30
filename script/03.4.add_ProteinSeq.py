import pandas as pd
from Bio import SeqIO
from pathlib import Path

# Define input and output paths
fasta_file = "/home/cx264/rds/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/05.Jinpeng/04.AlphaFold3/01.input/02.faa"
input_csv = "alphafold3_summary.threshold_0.6.csv"
output_csv = "extended_alphafold3_summary.threshold_0.6.csv"

# Read FASTA file and store protein names and sequences
protein_data = []
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        protein_data.append({
            "name": record.id,
            "sequence": str(record.seq)
        })

# Read the existing summary CSV
df = pd.read_csv(input_csv)

# Initialize new columns
df["PlantProteinName"] = "N/A"
df["PlantProteinSequence"] = "N/A"
df["BacteriaProteinName"] = "N/A"
df["BacteriaProteinSequence"] = "N/A"

# Map protein names and sequences based on PlantProteinID and BacteriaProteinID
num_plant_proteins = 333  # Number of plant proteins in 02.faa
for idx, row in df.iterrows():
    try:
        # PlantProteinID is 1-based index in FASTA (1 to 333)
        plant_idx = int(row["PlantProteinID"]) - 1
        if 0 <= plant_idx < num_plant_proteins:
            df.at[idx, "PlantProteinName"] = protein_data[plant_idx]["name"]
            df.at[idx, "PlantProteinSequence"] = protein_data[plant_idx]["sequence"]
        else:
            print(f"Warning: PlantProteinID {row['PlantProteinID']} out of range")

        # BacteriaProteinID is offset by num_plant_proteins (334th is BacteriaProteinID 001)
        bact_idx = num_plant_proteins + int(row["BacteriaProteinID"]) - 1
        if 0 <= bact_idx < len(protein_data):
            df.at[idx, "BacteriaProteinName"] = protein_data[bact_idx]["name"]
            df.at[idx, "BacteriaProteinSequence"] = protein_data[bact_idx]["sequence"]
        else:
            print(f"Warning: BacteriaProteinID {row['BacteriaProteinID']} out of range")
    except (ValueError, IndexError) as e:
        print(f"Error processing row {idx}: {e}")

# Reorder columns
columns = [
    "PlantProteinID", "PlantProteinName", "PlantProteinSequence", "PlantProteinAlignInfo",
    "BacteriaProteinID", "BacteriaProteinName", "BacteriaProteinSequence", "BacteriaProteinAlignInfo",
    "iPTM", "pTM", "ranking_score"
]
df = df[columns]

# Save to new CSV
df.to_csv(output_csv, index=False)
print(f"Extended summary table saved to {output_csv}")
