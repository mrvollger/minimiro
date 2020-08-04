#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd 



def get_color(x):
	if( x >= .99):
		return("255,103,0")
	elif(x >= .98):
		return("204,204,0")
	elif(x >= 0.9):
		mingray = 105; maxgray = 220
		gray = int( mingray + (maxgray-mingray) * (  1  - (x-.9)*10  ) )
		return(f"{gray},{gray},{gray}")
	else:
		# return purple
		return("147,112,219")


SEDEF_HEADER = "chr1   start1  end1    chr2    start2  end2    name    score   strand1 strand2 max_len aln_len comment aln_len.1 indel_a indel_b alnB    matchB  mismatchB   transitionsB     transversions   fracMatch       fracMatchIndel  jck     k2K     aln_gaps        uppercaseA      uppercaseB      uppercaseMatches        aln_matches  aln_mismatches  aln_gaps.1        aln_gap_bases   cigar   filter_score".strip().split()
SEDEF_HEADER = "chr1   start1  end1    chr2    start2  end2    name    score   strand1 strand2 max_len aln_len comment indel_a indel_b alnB    matchB  mismatchB   transitionsB     transversions   fracMatch       fracMatchIndel  jck     k2K     aln_gaps        uppercaseA      uppercaseB      uppercaseMatches        aln_matches  aln_mismatches  aln_gaps.1        aln_gap_bases   cigar   filter_score    count_ovls      sat_bases       total_bases     sat_coverage".strip().split()
DROP = ["aln_len.1", "aln_gaps.1", "cigar", "comment"]
DROP = ["aln_gaps.1", "cigar", "comment",     "count_ovls", "total_bases", "sat_coverage"]

# global var for inputs
args=None 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help="positional input")
    parser.add_argument("output", help="positional input")
    parser.add_argument("filt", help="positional input")
    parser.add_argument("-s", "--symetric", help="make sedef output symetric, set to 0 to disable.", default = 1)
    parser.add_argument("-l", "--minlength", help=" ", type=int, default=1000)
    parser.add_argument("-i", "--minidentity", help=" ", type=float, default=0.9)
    parser.add_argument("--sat", help="Remove dups that are this fraction of sat or more", type=float, default=0.70)
    parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()

    df =pd.read_csv(args.infile, sep="\t", names=SEDEF_HEADER, header=None, comment="#")
    # filter out high sat regions
    df = df[df.sat_coverage <= args.sat]

    df.drop(DROP, axis=1, inplace=True)
    #print(df.shape, sum( df["aln_gaps.1"] == df["aln_gaps"]), sum(df.aln_len == df["aln_len"]) )
    df["color"] = df.fracMatch.map(get_color)

    # make the symetric ones
    if(args.symetric):
            df2 = df.copy()
            df2[["chr1", "start1", "end1"]] = df[["chr2", "start2", "end2"]]
            df2[["chr2", "start2", "end2"]] = df[["chr1", "start1", "end1"]]
            df2["strand2"] = df["strand1"]
            df2["strand1"] = df["strand2"]
            df = pd.concat([df,df2], ignore_index=True)

    # make bed 9 format for the browser 
    df.sort_values(by=["chr1", "start1"], inplace=True)
    bed9 = ["chr1", "start1", "end1", "name", "fakeScore", "strand1", "start1", "end1", "color"]
    df["name"] = df.chr2 + ":" + df.start2.astype(str) + "-" + df.end2.astype(str)
    df["fakeScore"] = 0
    extra = [ col for col in SEDEF_HEADER if col not in bed9 and col not in DROP]
    df = df[bed9 + extra]
    df.rename(columns={"chr1":"#chr1"}, inplace=True)


    sd = df.loc[(df.aln_len >= args.minlength) & (df.fracMatch >= args.minidentity)] 
    filt = df.loc[(df.aln_len < args.minlength) | (df.fracMatch < args.minidentity)] 
    sd.to_csv(args.output, index=False, sep="\t")
    filt.to_csv(args.filt, index=False, sep="\t")


"""
chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100. Read more here.
The 9 additional optional BED fields are:

	name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
	score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
	shade									 
	score in range		≤ 166	167-277	278-388	389-499	500-611	612-722	723-833	834-944	≥ 945
	strand - Defines the strand. Either "." (=no strand) or "+" or "-".
	thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
	thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
	itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
"""
