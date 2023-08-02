#!/usr/bin/env python

import argparse
import os
import subprocess

def handle_args():
    parser = argparse.ArgumentParser(description='Run BIONJ on a dataset')
    parser.add_argument('--vcf', type=str, help='VCF file to run BIONJ on',
                        required=True)
    parser.add_argument('--outprefix', type=str, help='Prefix for output files',
                        required=True)
    parser.add_argument('--r-script', type=str, help='Path to r script file (not Rscript executable)',
                        required=True, default="run_bionj.R")
    parser.add_argument('--plink', type=str, help='Path to PLINK executable',
                        default='plink')
    parser.add_argument('--cleanup', action='store_true', help='Remove intermediate files')
    return parser.parse_args()

def check_input_file(file):
    """Check that the input file exists"""
    try:
        assert os.path.exists(file), f'File does not exist: {file}'
    except AssertionError as e:
        print(e)
        return False
    return True

def check_executable_runs(cmd):
    """Check that "cmd" is installed and working"""
    try:
        subprocess.run([cmd],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print(f'{cmd} is not available')
        return False
    return True

def check_args(args):
    """Check that the arguments are valid"""
    checks = [
        check_input_file(args.vcf),
        check_executable_runs(args.r_script),
        check_executable_runs(args.plink),
        ]
    return all(checks)

def run_plink_distance_matrix(vcf, outprefix):
    result = subprocess.run(['plink',
                             '--vcf', vcf,
                             '--distance', '1-ibs', 'square',
                             '--chr-set',  '38',
                             '--real-ref-alleles',
                             '--keep-allele-order',
                             '--const-fid',
                             '--out', outprefix])
    return result.returncode == 0

def run_bionj(r_script, outprefix):
    check_input_file(f'{outprefix}.mdist')
    check_input_file(f'{outprefix}.mdist.id')
    result = subprocess.run([r_script,
                             '-d', f'{outprefix}.mdist',
                             '-i', f'{outprefix}.mdist.id',
                             '-p', outprefix])
    return result.returncode == 0

def main():
    args = handle_args()
    if not check_args(args):
        return 1
    if not run_plink_distance_matrix(args.vcf, args.outprefix):
        return 2
    if not run_bionj(args.r_script, args.outprefix):
        return 4
    if args.cleanup:
        for f in os.listdir():
            if (os.path.isfile(f)
                and f.startswith(args.outprefix)
                and not f.endswith('.evec')
                and not f.endswith('.eigenvec')
                and not f.endswith('.mdist')
                and not f.endswith('.mdist.id')
                and not f.endswith('_bionj.nwk')
                and not f.endswith('_bionj.pdf')
                and not 'vcf' in f):
                os.remove(f)

    return 0

if __name__ == '__main__':
    main()
