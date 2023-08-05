#!/usr/bin/env python
import argparse
import os
import subprocess
import tempfile

eigenparams_template = \
'''genotypename:	{PREFIX}.ped
snpname:	{PREFIX}.map
indivname:	{PREFIX}.ped
outputformat:	EIGENSTRAT
genotypeoutname:	{PREFIX}.geno
snpoutname:	{PREFIX}.snp
indivoutname:	{PREFIX}.ind
familynames:	NO
'''

pcaparams_template = \
'''genotypename:	{PREFIX}.geno
snpname:	{PREFIX}.snp
indivname:	{PREFIX}.ind
evecoutname:	{PREFIX}.pca.evec
evaloutname:	{PREFIX}.eval
altnormstyle:	NO
numoutevec:	10
numoutlieriter:	0
numoutlierevec:	3
outliersigmathresh:	6
qtmode:	0
numthreads:	{THREADS}
'''


def handle_args():
    parser = argparse.ArgumentParser(description='Run PCA on a dataset')
    parser.add_argument('--vcf', type=str, help='VCF file to run PCA on',
                        required=True)
    parser.add_argument('--outprefix', type=str, help='Prefix for output files',
                        required=True)
    parser.add_argument('--plink', type=str, help='Path to PLINK executable',
                        default='plink')
    parser.add_argument('--smartpca', type=str, help='Path to smartpca executable',
                        default='smartpca')
    parser.add_argument('--convertf', type=str, help='Path to convertf executable',
                        default='convertf')
    parser.add_argument('--cleanup', action='store_true', help='Remove intermediate files')
    parser.add_argument('--threads', type=int, help='Number of threads to use',
                        default=1)
    args = parser.parse_args()
    return args

def check_input_file(vcf):
    """Check that the input file exists"""
    try:
        assert os.path.exists(vcf), 'VCF file does not exist'
    except AssertionError as e:
        print(e)
        return False
    return True

def check_args(args):
    """Check that the arguments are valid"""
    checks = [
        args.threads > 0,
        check_input_file(args.vcf),
        check_executable_runs(args.plink),
        check_executable_runs(args.smartpca),
        check_executable_runs(args.convertf)
        ]
    return all(checks)

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

def convert_vcf_to_plink(vcf, outprefix):
    """Convert a VCF file to PLINK format"""

    result = subprocess.run(['plink',
                             '--vcf', vcf,
                             '--recode', 'ped',
                             '--const-fid',
                             '--chr-set', '38',
                             '--real-ref-alleles',
                             '--keep-allele-order',
                             '--out', outprefix])

    # PLINK uses -9 for missing genotypes, but convertf expects 3
    try:
        with open(f'{outprefix}.ped.tmp', 'w') as f:
            subprocess.run(['sed', 's/-9/3/', f'{outprefix}.ped'],
                           stdout=f)
        os.rename(f'{outprefix}.ped.tmp', f'{outprefix}.ped')
    except FileNotFoundError:
        print('Could not find PLINK output file')
        return False
    return result.returncode == 0

def convert_plink_to_eigenstrat(outprefix):
    """Convert a PLINK file to eigenstrat format"""
    with tempfile.NamedTemporaryFile(mode='w') as tmp:
        tmp.write(eigenparams_template.format(PREFIX=outprefix))
        tmp.flush()
        result = subprocess.run(['convertf', '-p', tmp.name])
    return result.returncode == 0

def run_smartpca(outprefix, threads):
    """Run smartpca"""
    with tempfile.NamedTemporaryFile(mode='w') as tmp:
        tmp.write(pcaparams_template.format(PREFIX=outprefix, THREADS=threads))
        tmp.flush()
        result = subprocess.run(['smartpca', '-p', tmp.name])
    return result.returncode == 0

def run_plink_pca(vcf, outprefix):
    """Run PCA using PLINK"""
    result = subprocess.run(['plink',
                             '--vcf', vcf,
                             '--pca',
                             '--const-fid',
                             '--chr-set', '38',
                             '--real-ref-alleles',
                             '--keep-allele-order',
                             '--out', outprefix])
    try:
        with open(f'{outprefix}.eigenvec.tmp', 'w') as f:
            subprocess.run(['sed', 's/^0 //', f'{outprefix}.eigenvec'],
                           stdout=f)
        os.rename(f'{outprefix}.eigenvec.tmp', f'{outprefix}.eigenvec')
    except FileNotFoundError:
        print('Could not find PLINK output file')
        return False

    return result.returncode == 0

def main():
    args = handle_args()
    if not check_args(args):
        return 1
    if not convert_vcf_to_plink(args.vcf, args.outprefix):
        return 2
    if not convert_plink_to_eigenstrat(args.outprefix):
        return 4
    if not run_smartpca(args.outprefix, args.threads):
        return 8

    if args.cleanup:
        for f in os.listdir():
            if (os.path.isfile(f)
                and f.startswith(args.outprefix)
                and not f.endswith('.evec')
                and not f.endswith('.eval')
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

