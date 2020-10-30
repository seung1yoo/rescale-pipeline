#!/usr/bin/env python

import sys
import os
import argparse
import string
import vcf

def variant_filter(input_file, output_file):
    vcf_reader = vcf.Reader(open(input_file, 'r'))
    vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)
    for record in vcf_reader:
        chrID = record.CHROM
        pos = record.POS
        rsID = record.ID
        ref = record.REF
        alt = record.ALT[0]
        qual = record.QUAL
        list_varfilter = record.FILTER
        varfilter = ','.join(list_varfilter)
#         print(varfilter)
        if varfilter != "":
            continue
#         if "germline" in varfilter:
#             continue
#         if "panel_of_normals" in varfilter:
#             continue
#         if "strand_bias" in varfilter:
#             continue

        for sample in record.samples:
            try:
                af = float(sample['AF'])
            except:
                af = float(sample['AF'][0])
        
        if af < 0.01:
            continue

        # write output
        vcf_writer.write_record(record)
    return

def usage():
    message='''
python %s

-i, --vcf : vcf file
-o, --output : default(output.xlsx)
''' %sys.argv[0]
    print(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--vcf')
    parser.add_argument('-o', '--output', default="output.xlsx")
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.vcf) > 0

    except:
        usage()
        sys.exit(2)

    variant_filter(args.vcf,  args.output)

if __name__ == '__main__':
    main()
