import csv
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, mode='r', newline='', encoding='utf-8') as csvfile, \
     open(output_file, mode='w', newline='', encoding='utf-8') as tsvfile:
    
    reader = csv.reader(csvfile)
    writer = csv.writer(tsvfile, delimiter='\t')

    for row in reader:
        writer.writerow(row)

print(f"Conversion complete. TSV file saved as: {output_file}")
