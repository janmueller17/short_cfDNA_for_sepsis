#!/usr/bin/env python3

# This script extends the regions in a BED file to a minimum length to enable usage with AME.

import sys

def extend_regions(input_bed, output_bed, min_length):
    with open(input_bed, 'r') as input_file, open(output_bed, 'w') as output_file:
        for line in input_file:
            if line.startswith('#'):
                output_file.write(line)
                continue
            
            fields = line.strip().split('\t')
            chrom, start, end = fields[:3]
            length = int(end) - int(start)
            if length < min_length:
                extension = ((min_length - length)+1) // 2
                start = int(start) - extension
                end = int(end) + extension
                
            output_file.write('\t'.join([chrom, str(start), str(end)] + fields[3:]) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extend_bed.py input.bed output.bed min_length")
        sys.exit(1)

    input_bed_file = sys.argv[1]
    output_bed_file = sys.argv[2]
    min_length = int(sys.argv[3])
    extend_regions(input_bed_file, output_bed_file, min_length)
