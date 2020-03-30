#!/usr/bin/env python
import argparse
import sys
import pandas as pd 

# global var for inputs
args=None 


def read_genes(f):
	df = pd.read_csv(f, sep = "\t")
	df['exonStarts'] = df['exonStarts'].str.strip(",").str.split(',')
	df['exonEnds'] = df['exonEnds'].str.strip(",").str.split(',')
	print(df[["exonStarts", "exonEnds"]])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional input")
	parser.add_argument("-s", "--string", help="string option")
	parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
	parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
	parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()
	read_genes(args.infile)
