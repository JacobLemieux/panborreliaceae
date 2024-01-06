### translates a pan genome
# lemieux@broadinstitute.org
# Dec 26 2023

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import random

parser = argparse.ArgumentParser(description = 'translate list of sequences in a multifasta file')
parser.add_argument('-i', '--infile', help = 'input filename (fasta)', required=True)
parser.add_argument('-o', '--outfile', help = 'output filename (fasta)', required=True)
args = vars(parser.parse_args())

seqs = []
counter = 0
for seq_record in SeqIO.parse(args['infile'], "fasta"):
	#print(seq_record.seq)
	#if counter < 10:
	#	print("ID: " + seq_record.id + "Description: " + seq_record.description)
	#	counter += 1
	seq_record.seq = seq_record.seq.translate()
	seq_record.id = ""
	seq_record.description = seq_record.description.split(" ")[0] + "|" + seq_record.description.split(" ")[1]
	#print(seq_record.seq)
	seqs.append(seq_record)
    #seqs.append(seq_record.description.split()[1])
SeqIO.write(seqs, args['outfile'], "fasta")
