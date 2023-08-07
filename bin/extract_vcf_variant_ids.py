#!/usr/bin/env python3

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output IDS file (txt)')
    args = parser.parse_args()
    return args

def write_id_fields(vcf_in, ids_out):
    with open(ids_out, 'w') as f:
        result = subprocess.run(['bcftools', 'query',
                        '-f', '%ID\n',
                        vcf_in], stdout=f, stderr=subprocess.DEVNULL)
    return result.returncode

def main():
    args = cli()
    if not os.path.exists(args.input):
        print('Input file does not exist!')
        sys.exit(1)
    retcode = write_id_fields(args.input, args.output)
    return retcode


if __name__ == '__main__':
    import argparse
    import os
    import subprocess
    import sys
    main()
