import os
import json
import pandas as pd
import re
from pathlib import Path

# Define input and output paths
interaction_list = "002.interaction.threshold_0.6.list"
output_dir = "./03.output"
output_csv = "alphafold3_summary.threshold_0.6.csv"

# Function to extract first UniRef annotation from unpairedMsa
def extract_first_uniref(msa_content):
    # Split MSA into lines and find the first line starting with ">UniRef"
    lines = msa_content.split("\n")
    for line in lines:
        if line.startswith(">UniRef"):
            return line.strip()
    return "N/A"

# Read interaction list
with open(interaction_list, "r") as f:
    interactions = [line.strip() for line in f if line.strip()]

# Initialize list to store results
results = []

# Process each interaction
for pair in interactions:
    plant_id, bact_id = pair.split("_")
    pair_dir = os.path.join(output_dir, plant_id, f"{plant_id}_{bact_id}")
    
    # Load summary_confidences.json
    summary_file = os.path.join(pair_dir, f"{plant_id}_{bact_id}_summary_confidences.json")
    try:
        with open(summary_file, "r") as f:
            summary = json.load(f)
        iptm = summary.get("iptm", "N/A")
        ptm = summary.get("ptm", "N/A")
        ranking_score = summary.get("ranking_score", "N/A")
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error reading {summary_file}: {e}")
        iptm, ptm, ranking_score = "N/A", "N/A", "N/A"

    # Load data.json
    data_file = os.path.join(pair_dir, f"{plant_id}_{bact_id}_data.json")
    try:
        with open(data_file, "r") as f:
            data = json.load(f)
        # Extract unpairedMsa for plant (sequences[0]) and bacteria (sequences[1])
        plant_msa = data.get("sequences", [{}])[0].get("protein", {}).get("unpairedMsa", "")
        bact_msa = data.get("sequences", [{}])[1].get("protein", {}).get("unpairedMsa", "")
        plant_align_info = extract_first_uniref(plant_msa)
        bact_align_info = extract_first_uniref(bact_msa)
    except (FileNotFoundError, json.JSONDecodeError, IndexError) as e:
        print(f"Error reading {data_file}: {e}")
        plant_align_info, bact_align_info = "N/A", "N/A"

    # Append to results
    results.append({
        "PlantProteinID": plant_id,
        "PlantProteinAlignInfo": plant_align_info,
        "BacteriaProteinID": bact_id,
        "iPTM": iptm,
        "pTM": ptm,
        "ranking_score": ranking_score,
        "BacteriaProteinAlignInfo": bact_align_info
    })

# Create DataFrame and save to CSV
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)
print(f"Summary table saved to {output_csv}")
