#!/usr/bin/env python3

def cli():
    """
    Command line interface
    Usage:
        make_bootstrap_vcf.py --input-ids[FILE] --input-vcf[FILE] --output-vcf[FILE] --seed[INT]
        where input_ids is a file with one sample ID per line,
        input_vcf is a VCF file, seed is a random seed, and output_vcf is the output VCF file.
        Things will go wrong if input_ids is not derived from the input_vcf file, e.g. by running extract_sample_ids.py,
        but this is not checked.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-ids', required=True)
    parser.add_argument('--input-vcf', required=True)
    parser.add_argument('--output-vcf', required=True)
    parser.add_argument('--seed', required=True, type=int)
    args = parser.parse_args()
    return args

def subsample_ids(input_ids):
    with open(input_ids) as f:
        ids = [line.rstrip() for line in f]
    return Counter(random.choices(ids , k=len(ids)))

def main():
    args = cli()

    if not os.path.isfile(args.input_ids):
        sys.exit('Error: input_ids not found: ' + args.input_ids)
    if not os.path.isfile(args.input_vcf):
        sys.exit('Error: input_vcf not found: ' + args.input_vcf)
    if os.path.isfile(args.output_vcf):
        sys.exit('Error: output_vcf already exists: ' + args.output_vcf)

    random.seed(args.seed)
    ids = subsample_ids(args.input_ids)
    written = 1
    with tempfile.NamedTemporaryFile() as outf:
        with gzip.open(args.input_vcf) as f:
            for line in f:
                if line.startswith(b'#'):
                    outf.write(line)
                else:
                    line = line.decode()
                    fields = line.split('\t')
                    if fields[2] in ids:
                        for _ in range(ids[fields[2]]):
                            newline = '\t'.join(fields[:2] + [f'{fields[2]}.{written}'] + fields[3:])
                            written += 1
                            outf.write(newline.encode())
        outf.flush()
        subprocess.run(['bcftools', 'sort', '-Oz', '-o', args.output_vcf, outf.name])
        subprocess.run(['bcftools', 'index', '-f', args.output_vcf])

if __name__ == '__main__':
    import argparse
    import gzip
    import os
    import random
    import subprocess
    import sys
    import tempfile
    from collections import Counter
    main()
