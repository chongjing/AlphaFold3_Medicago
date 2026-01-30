###This script is to extract additional details like "unpairedMsa", "pairedMsa", and "templates" from AlphaFold3 pipeline processed individual protein
### and then merged to JSON file, for AlphaFold3 inference
### additional1.json and additional2.json are first and second protein, a pair for predicting interaction
### Chongjing Xia, xiachongjing@gmail.com

import json
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 code.py base.json additional1.json additional2.json ...")
        return

    base_file = sys.argv[1]
    
    # Load base JSON
    with open(base_file, 'r') as f:
        base_data = json.load(f)
    
    base_proteins = base_data.get('sequences', [])
    total_base = len(base_proteins)
    
    # Process additional files in sequence order
    for idx, add_file in enumerate(sys.argv[2:]):
        if idx >= total_base:
            break  # Ignore extra files if base proteins exhausted
            
        with open(add_file, 'r') as f:
            add_data = json.load(f)
            
            # Extract first protein from additional file
            try:
                add_protein = add_data['sequences'][0]['protein']
            except (KeyError, IndexError):
                continue
            
            # Get corresponding base protein
            target_protein = base_proteins[idx]['protein']
            
            # Merge template arrays
            if 'templates' in add_protein:
                target_protein.setdefault('templates', []).extend(add_protein['templates'])
            
            # Update MSA fields
            for field in ['unpairedMsa', 'pairedMsa']:
                if field in add_protein:
                    target_protein[field] = add_protein[field]

    # Save modified base JSON
    with open(base_file, 'w') as f:
        json.dump(base_data, f, indent=2)
    print(f"Merged {len(sys.argv)-2} additional files into {base_file}")

if __name__ == "__main__":
    main()
