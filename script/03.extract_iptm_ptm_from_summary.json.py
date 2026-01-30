import os
import json
import csv
import glob
import sys

def main(output_file):
    # Find all matching JSON files
    json_files = glob.glob('02.output/*/*_summary_confidences.json')

    # Prepare data storage
    data_rows = []

    for file_path in json_files:
        # Extract ID from filename
        filename = os.path.basename(file_path)
        file_id = filename.split('_summary_confidences.json')[0]

        # Load JSON data
        with open(file_path, 'r') as f:
            try:
                json_data = json.load(f)

                # Extract required values
                data_rows.append((
                    file_id,
                    json_data['iptm'],
                    json_data['ptm'],
                    json_data['ranking_score']
                ))
            except (KeyError, json.JSONDecodeError) as e:
                print(f"Error processing {file_path}: {str(e)}")
                continue

    # Write to CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'iptm', 'ptm', 'ranking_score'])
        writer.writerows(data_rows)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 code.py output.csv")
        sys.exit(1)
    main(sys.argv[1])
