import json
from pathlib import Path

def find_high_score_files(base_dir="02.output"):
    matching_files = []
    
    # Recursively find all target JSON files
    for json_path in Path(base_dir).rglob("*_summary_confidences.json"):
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
                
                # Check if all required keys exist
                required_keys = {"iptm", "ptm", "ranking_score"}
                if not required_keys.issubset(data.keys()):
                    continue
                
                # Validate numerical values
                if all(float(data[key]) >= 0.8 for key in required_keys):
                    matching_files.append(str(json_path))
                    
        except (json.JSONDecodeError, KeyError, ValueError) as e:
            print(f"Skipping {json_path} due to error: {e}")
    
    return matching_files

# Usage
result = find_high_score_files()
print("Files meeting criteria:")
for file in result:
    print(file)
