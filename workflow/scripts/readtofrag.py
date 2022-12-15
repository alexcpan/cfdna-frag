import pybedtools as bt
import argparse
from multiprocessing import Pool
import glob

parser = argparse.ArgumentParser()
parser.add_argument('--bed_dir',help='tmp input dir bedfile',required=True)
parser.add_argument('--fasta',help='fasta',required=True)
args = parser.parse_args()

BED_DIR = args.bed_dir
FASTA   = args.fasta
print("BED_DIR",BED_DIR)
print("FASTA",FASTA)

def get_nuc(bed):
    a = bt.BedTool(bed)\
        .nucleotide_content(fi=FASTA)\
        .saveas(bed+".nc")
        # .nucleotide_content(fi="/data/panalc/cfdna-frag/MPNST_frag/inputs/hg38.fa")\

files = glob.glob(BED_DIR + "/*.pe.*")
print("Files to merge: ", len(files))

pool = Pool() #note the default will use the optimal number of workers
pool.map(get_nuc,files)

print("Finished merging")
