from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import random
import pandas as pd

parser = argparse.ArgumentParser(description = 'translate list of sequences in a multifasta file')
parser.add_argument('-i', '--infile', help = 'input filename (SpLip output)', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename (tsv)', required=True)
args = vars(parser.parse_args())

lipo_dict = {}

splip_file = open(args['infile'], 'r')
lines = splip_file.readlines()
for line in lines:
	if line[0] == ">":
		line = line.strip()
		lipo_dict[line.split(":")[0].split(" ")[-1].split("|")[1]] = line.split(":")[1]

# write this out to a csv
df = pd.DataFrame.from_dict(lipo_dict, orient = 'index')
df.to_csv(args['outfile'])