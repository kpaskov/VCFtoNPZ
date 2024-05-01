import numpy as np
from scipy.sparse import csc_matrix, save_npz, hstack
import time
import argparse
import gzip
from pysam import VariantFile, TabixFile
import json
import os
import itertools

t0 = time.time()

parser = argparse.ArgumentParser(description='Pull genotypes.')
parser.add_argument('vcf_file', type=str, help='VCF file to pull from.')
parser.add_argument('assembly', type=str, help='Human genome reference used.')
parser.add_argument('data_dir', type=str, help='Data directory.')
parser.add_argument('name', type=str, help='Name of filter.')
parser.add_argument('chrom', type=str, help='Chromosome of interest.')
parser.add_argument('--DP', type=int, default=None, help='DP (depth) cutoff for filter.')
parser.add_argument('--GQ', type=int, default=None, help='GQ (genotype quality) cutoff for filter.')
parser.add_argument('--AB', type=flaot, default=None, help='AB (allele balance) cutoff for filter.')
parser.add_argument('--maxsize', type=int, default=10000000000, help='Amount of memory per block.')
parser.add_argument('--id_mapper_file', type=str, default=None, help='File that maps old ids to new ones.')
parser.add_argument('--id_mapper_sep', type=str, default='\t', help='Separater to parse id_mapper_file.')
parser.add_argument('--old_id_index', type=int, default=0, help='Index of old_id in id_mapper_file.')
parser.add_argument('--new_id_index', type=int, default=1, help='Index of new_id in id_mapper_file.')
args = parser.parse_args()

# depth of R 10x, a genotype quality of R 25, and a ratio of alternative allele reads/total reads R 0.2. 

if not os.path.exists('%s/filters/%s' % (args.data_dir, args.name)):
    os.makedirs('%s/filters/%s' % (args.data_dir, args.name))

with open('%s/genotypes/info.json' % args.data_dir, 'r') as f:
    info = json.load(f)

chrom_int = 23 if args.chrom == 'X' else 24 if args.chrom == 'Y' else 25 if args.chrom == 'MT' else int(args.chrom)

gen_mapping = {'./.': -1, '0/0': 0, '0|0': 0, '0/1': 1, '0|1': 1, '1/0': 1, '1|0': 1, '1/1': 2, '1|1': 2}

def process_header(vcf):
    sample_ids = [x.replace('.', '_') for x in vcf.header.samples]

    if args.id_mapper_file is not None:
        old_id_to_new_id = dict()
        with open(args.id_mapper_file, 'r') as f:
            for line in f:
                pieces = line.strip().split(args.id_mapper_sep)
                if len(pieces)>args.old_id_index and len(pieces)>args.new_id_index:
                    old_id_to_new_id[pieces[args.old_id_index]] = pieces[args.new_id_index]
        sample_ids = [old_id_to_new_id[x] for x in sample_ids]

    sample_file = '%s/genotypes/samples.json' % args.data_dir
    with open(sample_file, 'r') as f:
        stored_sample_ids = json.load(f)
        assert sample_ids == stored_sample_ids

    return sample_ids, vcf.header.contigs


def process_body(records, sample_ids):

    data, indices, indptr, index = np.zeros((args.maxsize,), dtype=bool), np.zeros((args.maxsize,), dtype=int), [0], 0

    for line in records:
        pieces = line.strip().split('\t')
        fmt = pieces[8].strip().split(':')

        # pull chrom_coord information
        pos, _, ref, alt = pieces[1:5]

        # pull genotypes
        gen_index = fmt.index('GT')
        dp_index = fmt.index('DP')
        ab_index = fmt.index('AB')
        gq_index = fmt.index('GQ')
        for i, piece in enumerate(pieces[9:]):
            segment = piece.split(':', maxsplit=gen_index+1)
            if segment[gen_index] == '0/1' or segment[gen_index] == '0|1' or segment[gen_index] == '1|0':
                is_ok = True
                if (args.DP is not None) and ((segment[dp_index] == '.') or (int(segment[dp_index]) < args.DP)):
                    is_ok = False
                elif (args.AB is not None) and ((segment[ab_index] == '.') or (float(segment[ab_index]) < args.AB)):
                    is_ok = False
                elif (args.GQ is not None) and ((segment[gq_index] == '.') or (int(segment[gq_index]) < args.GQ)):
                    is_ok = False
                
                if not is_ok:
                    indices[index] = i
                    data[index] = True
                    index += 1
        indptr.append(index)

    fil = csc_matrix((data[:index], indices[:index], indptr), shape=(len(sample_ids), len(indptr)-1), dtype=bool)

    # Save to file
    save_npz('%s/filters/%s/chr.%s.%d.gen' % (args.data_dir, args.chrom, info['batch_num']), fil)

with open('%s/filters/%s/info.json' % (args.data_dir, args.name), 'w+') as f:
    json.dump(vars(args), f, indent=4)

vcf = VariantFile(args.vcf_file)
sample_ids, contigs = process_header(vcf)

contig = None
if args.chrom in contigs:
    contig = contigs[args.chrom]
elif 'chr%s' % args.chrom in contigs:
    contig = contigs['chr%s' % args.chrom]
else:
    raise Exception('Trouble finding contig', args.chrom, 'in', contigs)
print('Chrom length', contig.length)


vcf_files = [args.vcf_file]

if np.all([os.path.isfile(vcf_file + '.tbi') for vcf_file in vcf_files]):
    vcfs = [TabixFile(vcf_file, parser=None) for vcf_file in vcf_files]

    if info['batch_size'] != -1:
        start_pos, end_pos = info['batch_num']*info['batch_size'], (info['batch_num']+1)*info['batch_size']
        print('Interval', start_pos, end_pos)
        if start_pos < contig.length:
            process_body(itertools.chain(*[vcf.fetch(reference=contig.name, start=start_pos, end=end_pos) for vcf in vcfs]), sample_ids)
        else:
            print('Interval (%d-%d) is longer than chromosome (length=%d).' % (start_pos, end_pos, contig.length))
    else:
        process_body(itertools.chain(*[vcf.fetch(reference=contig.name) for vcf in vcfs]), sample_ids)
else:
    print('Error, .tbi files are missing.')

print('Completed in ', time.time()-t0, 'sec')

