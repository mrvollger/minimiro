#!/usr/bin/env python
import argparse
import sys
import re 
import pandas as pd

rm_color_map = {
	"DNA": "Pink",
	"LINE": "DkGray",
	"Low_complexity": "Yellow",
	"LTR": "FGreen",
	"RC": "Green",
	"rRNA": "PwdBlue",
	"Satellite": "Red",
	"scRNA": "SkyBlue",
	"Simple_repeat": "Orange",
	"SINE": "CbBlue",
	"snRNA": "Blue",
	"srpRNA": "DkBlue",
	"tRNA": "Purple",
	"Unknown": "LtGray" }
COLORS=list(rm_color_map.values())

parser = argparse.ArgumentParser(description="")
parser.add_argument("-r",  help="input ref file (fasta/fastq)" )
parser.add_argument("-q",  help="input query file (fasta/fastq)", type=argparse.FileType('r'))
parser.add_argument("--paf", help="plot PAF alignment from minimap2, must have cs string")
parser.add_argument("-o", "--out", help="output file", type=argparse.FileType('w') )
parser.add_argument("--rm", nargs="+", help="repeatmasker .out file, can be a list of files", default=[])
parser.add_argument("--dm", nargs="+", help="dupmakser file, parsed to add color,  can be a list of files", default=[])
parser.add_argument("--bed", nargs="+", help="Must be bed 6 format", default=[])
parser.add_argument("-x", "--preset", help="[asm20]", default="map-ont")
parser.add_argument("-t", "--threads", help="[4]", type=int, default=4)
parser.add_argument("-n", "--bestn", help="[100]", type=int, default=100)
parser.add_argument("--start", help="start positioon", type=int, default=0)
parser.add_argument("-e", "--exclude", nargs="+", help="exclude these conitg names from anlysis", default=[])
parser.add_argument("-s", "--score", help="threashold for reporting alignments, same as -s in miropeats [100]", type=int, default=100)
parser.add_argument("-d", "--drop", help="if the alignemnt drops by [score/10] terminate the segment", type=int, default=None)
args = parser.parse_args()
if(args.drop is None):
	args.drop = int(args.score / 10)

# patter for cs string in an alignment
pattern=re.compile("(:)([0-9]+)|(\*)([a-z][a-z])|([=\+\-])([A-Za-z]+)")
consume_r = set([":", "*", "-"])
consume_q = set([":", "*", "+"])
minus = set(["*", "+", "-"])

# turn cs string inot lsit of tupples 
def parse_cs(cs):
	rtn = []
	for match in re.findall(pattern, cs):
		match = [x for x in match if x!= '']
		assert len(match) == 2, match 
		if(match[0] == ":"):
			match[1] = int(match[1])
		rtn.append( (match[0], match[1]) )
	return(rtn)

def high_scoring_segs(hit, cs, q_ctg, q_ctg_len):
	rtn = []

	r_st = hit.r_st
	q_st = hit.q_st; q_en = hit.q_en
	if(hit.strand == -1 ): q_st = hit.q_en; q_en = hit.q_st

	score = 0
	i_max = 0
	r_seg_start = r_st
	r_seg_end = r_st
	q_seg_start = q_st
	q_seg_end = q_st
	r_cur = r_st
	q_cur = q_st
	
	# if we are passing a paf, plot them exactly and do not apply the segment algorithum
	if(args.paf):
		rtn.append((hit.ctg, hit.r_st, hit.r_en, hit.ctg_len, q_ctg, hit.q_st, hit.q_en, q_ctg_len))
		return(rtn)

	total = len(cs)
	progress = 0
	for opt, val in cs:
		progress += 1
		sys.stderr.write("\rParssing {} to {}: {:.2%}".format(q_ctg, hit.ctg, progress/total ))
		# length of operation
		opt_l = 0
		if(opt == ":"):
			opt_l = val
		elif(opt == "*"):
			opt_l = 1
		elif(opt in ["-", "+"] ):
			opt_l = len(val)
		
		# iterate over operations
		for i in range(opt_l):
			# update score 
			if(opt == ":"):
				score += 1
			else: #elif( opt in minus ):
				score -= 1
			
			# incremen current d segment 
			if(score >= i_max):
				i_max = score
				r_seg_end = r_cur
				q_seg_end = q_cur		

			#print(score, i_max)
			#if(score <= 0 or (r_cur + 1) == hit.r_en  ):
			if(score < 0 or score <= (i_max - args.drop) or (r_cur + 1) == hit.r_en  ):
				#print(score, i_max)
				if( i_max >= args.score  ):
					rtn.append((hit.ctg, r_seg_start, r_seg_end, hit.ctg_len, q_ctg, q_seg_start, q_seg_end, q_ctg_len))
				i_max = 0
				score = 0
				r_seg_start = r_cur+1
				r_seg_end = r_cur+1
				q_seg_start = q_cur+1
				q_seg_end = q_cur+1
			
			# incrememnt positions:
			if(opt in consume_r):
				r_cur += 1
			if(opt in consume_q):
				q_cur += 1*hit.strand

	sys.stderr.write("\n")
	assert r_cur == hit.r_en
	assert q_cur == q_en 
	return(rtn)




"""
Col	Type	Description
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start coordinate (0-based)
4	int	Query end coordinate (0-based)
5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start coordinate on the original strand
9	int	Target end coordinate on the original strand
10	int	Number of matching bases in the mapping
11	int	Number bases, including gaps, in the mapping
12	int	Mapping quality (0-255 with 255 for missing)
"""
"""
An Alignment object can be converted to a string with str() in the following format:

	q_st  q_en  strand  ctg  ctg_len  r_st  r_en  mlen  blen  mapq  cg:Z:cigar_str
"""

segs = []

if(args.paf):
	#
	# READ IN MINIMAP2 PAF
	#
	firstcol = ["q_ctg", "q_len", "q_st", "q_en", "strand", "ctg", "ctg_len", "r_st", "r_en", "match", "span", "mapq"]
	pd_in = [[] for x in range(len(firstcol)+1)]
	for line in open(args.paf):
		tokens = line.strip().split()
		# skip alignments with low scores
		if(int(tokens[9]) < args.score ):
			continue 
		# parse the paf input 
		for idx, token in enumerate(tokens):
			if(idx < len(firstcol)):
				try:
					token = int(token)
				except ValueError:
					pass
				pd_in[idx].append(token)
			else:
				match = re.match("cs:Z:\S+", token)
				if(match):
					pd_in[len(firstcol)].append(token)
					break
	
	paf = pd.DataFrame(pd_in)
	paf = paf.transpose()
	paf.columns = firstcol + ["cs"]
	paf["strand"].replace({"+":1, "-":-1}, inplace=True)
	
	#print(paf) 
	for idx, row in paf.iterrows():
		cs = re.match("cs:Z:(\S+)", row["cs"])
		assert cs, "Must have cs string, line {}".format(idx+1)
		cs = parse_cs( cs.groups()[0] )
		segs += high_scoring_segs(row, cs, row["q_ctg"], row["q_len"] )
		#segs += high_scoring_segs(hit, cs, name, qlen )

elif(args.r):
	#
	# RUN MINIMAP
	#
	import mappy as mp
	ref = mp.Aligner(args.r, preset = args.preset, best_n=args.bestn )  # load or build index
	if not ref: raise Exception("ERROR: failed to load/build index")

	sys.stderr.write("\rLoaded index\n")

	# set minimap2 alignemtn flag
	allchain = 0x800000
	nodiag = 0x001
	flag = allchain + nodiag

	# run alignemtns across inputs
	counter = 0
	for name, seq, qual in mp.fastx_read(args.q.name):
		qlen = len(seq)	
		if(name in args.exclude):
			continue 
		sys.stderr.write("\rAligning {}".format(name))
		for hit in ref.map(seq, cs=True, extra_flags=flag):
			counter += 1; sys.stderr.write("\rProccessed {} alignments".format(counter))
			cs = parse_cs(hit.cs)
			qlen 
			segs += high_scoring_segs(hit, cs, name, qlen )
			#for seg in segs: print(seg)
			#print(segs)
else:
	sys.stderr.write("Needs ref and query or paf input!")
	exit(1)

#
# PARSE MINIMAP FOR MIROPEATS
#
segs = pd.DataFrame(segs, columns=["r_ctg", "r_st", "r_en", "r_len", "q_ctg", "q_st", "q_en", "q_len"]).sort_values(by=["r_ctg", "r_st"])

# remove things beofre start 
segs.r_st = segs.r_st - args.start
segs.r_en = segs.r_en - args.start
segs.r_len = segs.r_len - args.start

segs.q_st = segs.q_st - args.start
segs.q_en = segs.q_en - args.start
segs.q_len = segs.q_len - args.start

segs = segs[ (segs.r_en >= 0) & (segs.q_en >= 0) ]
segs.loc[segs.r_st < 0, ["r_st"]] = 0 
segs.loc[segs.q_st < 0, ["q_st"]] = 0 

# set up contig sizes 
contigs = pd.DataFrame( [ tuple(x) + (1,) for x in segs[["r_ctg","r_len"]].drop_duplicates().values] + [ tuple(x) +(0,) for x in segs[["q_ctg","q_len"]].drop_duplicates().values], columns=["ctg", "len", "ref"] )
contigs.sort_values(by=["ref", "len"], inplace=True)
contigs.drop_duplicates(subset=["ctg"], inplace=True, keep="last")
def get_contig_num(contig):
	return(list(contigs["ctg"]).index(contig))

def get_contig_len(contig):
	return( contigs["len"][contigs["ctg"] == contig ].iloc[0] )

contigs["y"] = contigs.ctg.map(get_contig_num)

# add contig lines and labels 
CTGS = ""
for idx, row in contigs.iterrows():
	CTGS += "drop {} mul ({}) drop {} mul ({}) printnames\n".format(row.y, row.ctg, row.y, row.ctg)
	CTGS += "{} xstretch div cm  drop {} mul tagends\n".format(row["len"], row.y)
	if(len(args.rm)> 0):
		CTGS += "{} (RepeatMasker) printname_right\n".format(row.y + 0.02 + 0.035)
	if(len(args.dm)> 0):
		CTGS += "{} (DupMasker) printname_right\n".format(row.y + 0.08 + 0.035)

# add homology lines
segs.insert(loc=0, column='id1', value=segs.r_ctg.map(get_contig_num))
segs.insert(loc=5, column='id2', value=segs.q_ctg.map(get_contig_num))

segs["r_ctg"] = "(" + segs["r_ctg"] + ")"
segs["q_ctg"] = "(" + segs["q_ctg"] + ")"
segs["func"] = "printcontig"

SEGS=""
counter = 0
for idx, row in segs.iterrows():
	
	if(row.r_ctg == row.q_ctg and row.r_st == row.q_st and row.r_en == row.q_en):
		#print("skipping:{}".format(row, file=sys.stderr))
		continue

	SEGS += COLORS[counter] + "\n"
	for x in row:
		SEGS += str(x) + " "
	SEGS += "\n"
	
	# cycle the colors
	counter += 1
	if(counter == len(COLORS) ):
		counter = 0
SEGS += "Black\n"


#
# READ IN REPEATMASKER
#
def get_rm_end(strand):
	if(strand in ["+", "F"]):
		return("RIGHT_end")
	return("LEFT_end")
def get_rm_color(rclass):
	if(rclass in rm_color_map):
		return(rm_color_map[rclass])
	return("Black")

if(len(args.rm) > 0):
	rms = []
	for frm in args.rm:
		rms.append( pd.read_csv(frm, delim_whitespace=True, header=None, skiprows=[0,1,2], comment="*", 
			names = ["score", "div", "del", "ins", "q_ctg", "q_st", "q_en", "q_left", "strand", "repeat", "class", "r_st", "r_en", "r_left", "id"]) )
		
	rm = pd.concat(rms, ignore_index=True) 
	rm = rm[rm.q_ctg.isin(contigs.ctg)]

	rm.q_st = rm.q_st - args.start
	rm.q_en = rm.q_en - args.start
	rm = rm[rm.q_st >= 0]

	rm["y"] = rm.q_ctg.map(get_contig_num) + 0.02
	rm["func"] = rm.strand.map(get_rm_end)
	rm["rclass"] = rm["class"].str.split("/|\?", expand=True)[0]
	rm["color"] = rm.rclass.map(get_rm_color)

	RM = ""
	for idx, row in rm.iterrows():
		if(row.q_en < get_contig_len(row.q_ctg) ):
			RM += "{}\n{} {} {} {}\n".format(row.color, row.y, row.q_st, row.q_en, row.func)
	
	# add a color legend
	RM += "0 -1.5 cm moveto\nBlack\n(RepeatMasker:\t\t) 50 string cvs show\n"
	for rclass in rm_color_map:
		color = rm_color_map[rclass]
		RM += "{}\n({}\t\t) 50 string cvs show\n".format(color, rclass)
	RM += "Black\n"
else:
	RM = ""

#
# READ in DupMasker
#
def hex_to_rgb(h):
	rgb= tuple(int( h.lstrip("#")[i:i+2], 16)/255 for i in (0, 2, 4))
	rgb = "{} {} {} setrgbcolor\n".format(rgb[0], rgb[1], rgb[2])
	return(rgb)

if(len(args.dm) > 0):
	dms = []
	for fdm in args.dm:
		dms.append( pd.read_csv(fdm, delim_whitespace=True))
		
	dm = pd.concat(dms, ignore_index=True) 
	if("qChr" in dm):
		dm["chr"] = dm["qChr"]
		dm["chrStart"] = dm["qStart"]
		dm["chrEnd"] = dm["qEnd"]
		dm["orient"] = dm["Orient"]

	dm = dm[dm["chr"].isin(contigs.ctg)]
	
	dm.chrStart = dm.chrStart - args.start
	dm.chrEnd = dm.chrEnd - args.start
	dm = dm[dm.chrStart >= 0]

	dm["y"] = dm["chr"].map(get_contig_num) + 0.08
	dm["func"] = dm.orient.map(get_rm_end)

	DM = ""
	for idx, row in dm.iterrows():
		if(row.chrEnd < get_contig_len( row["chr"] ) ):
			DM += hex_to_rgb(row.color)
			DM += "{} {} {} {}\n".format(row.y, row.chrStart, row.chrEnd, row.func)
	DM += "Black\n"
else:
	DM = ""

#
# Read in genes
#
def get_bed_end(strand):
	if(strand in ["+", "F"]):
		return("RIGHT_end")
	return("LEFT_end")

if(len(args.bed) > 0):
	beds = []
	for fbed in args.bed:
		beds.append( pd.read_csv(fbed, delim_whitespace=True, names = ["chr", "start", "end", "name", "score", "strand"]))
	bed = pd.concat(beds, ignore_index=True) 
	bed = bed[bed["chr"].isin(contigs.ctg)]
	bed["y"] = bed["chr"].map(get_contig_num) + 0.14
	bed["func"] = bed.strand.map(get_bed_end)
	
	BED = ""
	offset = 0
	for name, group in bed.groupby("name"):
		BED += "{} {} {} {} draw_line\n".format(min(group.start), row.y+offset + 0.025 , max(group.end), row.y+offset+ 0.025)
		for idx, row in group.iterrows():
			if(row.end < get_contig_len( row["chr"] ) ):
				BED += "{} {} {} {}\n".format(row.y + offset , row.start, row.end, row.func)
		
		BED += "{} {} ({}) printname_right_2\n".format( max(group.end), row.y + offset + 0.015, name)
		
		offset += 0.03
	#print(bed)
	#print(bed[bed["name"] =="AK6"])
else:
	BED=""


#
# SET SOME PS CONSTANTS 
#
# Set dimensions in centimetres
cm_to_px = 28.3
rightmost = max(contigs["len"])
contig_count = len( contigs.index )
print("CONTIG COUNT:{}".format(contig_count, file=sys.stderr))
LEFT_MARGIN = 2 ; 
BOT_MARGIN = 2 ; 
GRAPH_MARGIN = 20 ;
NAME_MARGIN = 0.6 ;
PAGE_LENGTH = 7 * ( contig_count ) ;
# Time to play with some fixed point arithmetic
contig_stretch100 =int( ( rightmost * 100 )/ GRAPH_MARGIN );
repeat_count = len(segs.index)
NAME_FONT_SIZE = 8;
REP_WIDTH = 0.25;


# add a scale
SCALE_WIDTH = 0.1;                                                                                                                                                                                           
SCALE_BELOW = -0.3;
SCALE_NUM_BELOW = -.5 #str( SCALE_WIDTH * SCALE_BELOW ) ;
TICK_WIDTH = 0;
SCALE_FONT_TYPE = "Helvetica";
SCALE_FONT_SIZE = 7;

digits = len(str(rightmost))
MAJ_SCALE_STEP = 10**(digits - 1)
if(MAJ_SCALE_STEP < 100):
    MAJ_SCALE_STEP = 100;
MIN_SCALE_STEP = int( MAJ_SCALE_STEP / 10); 

HEIGHT = int(PAGE_LENGTH * cm_to_px*1.2)
WIDTH = int( (GRAPH_MARGIN+2*LEFT_MARGIN) * cm_to_px*1.2)  

# y postions from contig name
def y_pos(contig):
	#/drop {PAGE_LENGTH} {contig_count} div cm def
	#/y1 drop number1 mul def
	return( PAGE_LENGTH/contig_count * get_contig_num(contig) * cm_to_px ) 

# x postions from contig position
def x_pos( contig_pos ):
	#contig_stretch100 =int( ( rightmost * 100 )/ GRAPH_MARGIN );"
	#/xstretch {contig_stretch100} 100 div def"
    #/x1 start1 xstretch div cm def
	return(contig_pos / (rightmost / GRAPH_MARGIN) * cm_to_px)

def printcontig(contig1, start1, end1, contig2, start2, end2, color = None):
	rtn="" # postscript to be executed 	
	if(contig1 == contig2 and start1 == start2 and end1 == end2):
		return(rtn)
	x1 = x_pos(start1); x2 = x_pos(end1); x3 = x_pos(start2); x4 = x_pos(end2)
	y1 = y_pos(contig1); y2 = y_pos(contig2)
	
	middle_y1 = (y1+y2)/2
	if(y1 == y2):
		middle_y1 = middle_y1 + PAGE_LENGTH / contig_count * 0.75
	middle_y2 = middle_y1

out = f"""
%\!PS-Adobe-1.0
%%Creator: mvollger
%%Title: Repeats Graphic
%%CreationDate: DATE
%%Pages: 1
%%DocumentFonts: Helvetica-Bold Helvetica
%%BoundingBox: 0 0 612 792
%%EndComments
<< /PageSize [{WIDTH} {HEIGHT}] >> setpagedevice
/cm {{28.35 mul}} def
%%EndProlog
%%Page: 0 1

/repwidth {REP_WIDTH} cm def
/linkwidth 0.5 def
/gapwidth 0 def
/leftmargin {LEFT_MARGIN} cm def
/bottommargin {BOT_MARGIN} cm def
/pageheight {PAGE_LENGTH} cm def
leftmargin bottommargin translate
/Helvetica-Bold findfont 13 scalefont setfont
0 pageheight 0.0 cm add moveto
%(Threshold Score = {args.score}) show
/Helvetica findfont {NAME_FONT_SIZE} scalefont setfont


/printcontig {{
        /length2 exch def
        /end2 exch def
        /start2 exch def
        /name2 exch def
        /number2 exch def
        /length1 exch def
        /end1 exch def
        /start1 exch def
        /name1 exch def
        /number1 exch def
        /y1 drop number1 mul def
        /x1 start1 xstretch div cm def
        /x2 end1 xstretch div cm def
        /y2 drop number2 mul def
        /x3 start2 xstretch div cm def
        /x4 end2 xstretch div cm def
        /middley1 y1 y2 add 2 div def
        number1 number2 eq {{/middley1 middley1 drop 0.75 mul add def}} if
        /middley2 middley1 def
        number1 number2 eq {{start1 end1 sub start2 end2 sub ne start1 end2 ne and {{/middley1 middley1 drop 0.2 mul add def /middley2 middley2 drop 0.15 mul add def}} if }} if
        start1 start2 gt {{ /temp middley1 def /middley1 middley2 def /middley2 temp def }} if

        repwidth setlinewidth
        newpath
        x1 y1 moveto
        x2 y1 lineto
        stroke
        newpath
        x3 y2 moveto
        x4 y2 lineto
        stroke
        linkwidth setlinewidth
        newpath
        x1 y1 moveto
        x1 middley1 x3 middley1 x3 y2 curveto
        stroke
        newpath
        x2 y1 moveto
        x2 middley2 x4 middley2 x4 y2 curveto
        stroke
        %y1 name1 y2 name2 printnames
        %length1 xstretch div cm y1 tagends
        %length2 xstretch div cm y2 tagends
}} def

/graphicmargin {GRAPH_MARGIN} cm def
/printnames {{
        /name2 exch def
        /height2 exch def
        /name1 exch def
        /height1 exch def
        graphicmargin {NAME_MARGIN} cm add height1 moveto
        0 {NAME_FONT_SIZE} -2 div rmoveto
        name1 100 string cvs show
        graphicmargin {NAME_MARGIN} cm add height2 moveto
        0 {NAME_FONT_SIZE} -2 div rmoveto
        name2 100 string cvs show
        }} def

/printname_right {{
        /name1 exch def
        /height1 exch def
        /y1 drop height1 mul def
        graphicmargin 0.6 cm add y1 moveto
        0 8 -2 div rmoveto
        name1 100 string cvs show

}} def

/printname_right_2 {{
        /name1 exch def
        /height1 exch def
        /end1 exch def
        /y1 drop height1 mul def
        /x1 end1 xstretch div cm def
        0.1 x1 add y1 moveto
        %0 8 -2 div rmoveto
        name1 50 string cvs show

}} def


/printname_left {{
        /name1 exch def
        /height1 exch def
        /y1 drop height1 mul def
        0 y1 moveto
        name1 100 string cvs show

}} def

/draw_line {{                                                                                                                                               
        /height2 exch def
        /end2 exch def
        /height1 exch def
        /end1 exch def
        /y1 drop height1 mul def
        /x1 end1 xstretch div cm def
        /y2 drop height2 mul def
        /x2 end2 xstretch div cm def
        0.1 setlinewidth
        x1 y1 newpath moveto
        x2 y2 lineto
        stroke
}} def


/tagends {{
        /ypos exch def
        /xpos exch def
        newpath
        xpos ypos moveto
        xpos {REP_WIDTH} cm add ypos {REP_WIDTH} cm add lineto
        xpos {REP_WIDTH} cm add ypos {REP_WIDTH} cm sub lineto
        closepath
        fill
        newpath
        0 setlinewidth
        0 ypos moveto
        xpos ypos lineto stroke
        }} def

% ------------------------------------------------
% FUCNTIONS FROM THE ANNOTATED MIROPEATS PIPELINE

% Colors
/Pink {{0.85 0.45 0.6 setrgbcolor}}   def
/Red {{0.6 0 0 setrgbcolor}}        def
/Magenta {{1 .3 .5 setrgbcolor}}  def
/Orange {{1 .62 .14 setrgbcolor}} def
/Yellow {{0.9 0.8 0.05 setrgbcolor}}     def

/FGreen {{0.4 0.5 0.2 setrgbcolor }}    def
/Green {{0.0 0.6 0.0 setrgbcolor }}   def
/DrbGreen {{0.6 0.8 0.6 setrgbcolor }}    def
/PwdBlue {{0.7 0.7 1 setrgbcolor}}        def
/SkyBlue {{0.7 0.8 1 setrgbcolor}} def
/CbBlue {{0.6 0.6 1 setrgbcolor}} def
/Blue {{0 0.4 0.8 setrgbcolor}}   def
/DkBlue {{0 0 0.6 setrgbcolor }}  def
/Purple {{0.47 0 0.47 setrgbcolor }}      def

/LtGray {{.75 .75 .75 setrgbcolor}}       def
/LtGray2 {{.71 .71 .71 setrgbcolor}}       def
/MedGray {{.5 .5 .5 setrgbcolor}} def
/DkGray {{0.3 0.3 0.3 setrgbcolor}}       def
/White {{1 1 1 setrgbcolor }}     def
/Black {{0 0 0 setrgbcolor }}     def


/RIGHT_end
{{
   /end1 exch def
   /start1 exch def
   /number1 exch def
   /y1 drop number1 mul def
   /y2 y1 10 add def
   /x1 start1 xstretch div cm def
   /x2 end1 xstretch div cm def
   /ym y1 5 add def
   newpath x1 y1 moveto x2 ym lineto x1 y2 lineto  x1 y1 lineto
   closepath
   fill stroke

}} def

/RIGHT_end2
{{
   /end1 exch def
   /start1 exch def
   /number1 exch def
   /y1 drop number1 mul def
   /y2 y1 6 add def
   /x1 start1 xstretch div cm def
   /x2 end1 xstretch div cm def
   /ym y1 3 add def
   newpath x1 y1 moveto x2 ym lineto x1 y2 lineto  x1 y1 lineto
   closepath
   fill stroke

}} def

/LEFT_end
{{
   /end1 exch def
   /start1 exch def
   /number1 exch def
   /y1 drop number1 mul def
   /y2 y1 10 add def
   /x1 start1 xstretch div cm def
   /x2 end1 xstretch div cm def
   /ym y1 5 add def
   newpath x1 ym moveto x2 y2 lineto x2 y1 lineto  x1 ym lineto
   closepath
   fill stroke

}} def

/LEFT_end2
{{
   /end1 exch def
   /start1 exch def
   /number1 exch def
   /y1 drop number1 mul def
   /y2 y1 6 add def
   /x1 start1 xstretch div cm def
   /x2 end1 xstretch div cm def
   /ym y1 3 add def
   newpath x1 ym moveto x2 y2 lineto x2 y1 lineto  x1 ym lineto
   closepath
   fill stroke

}} def

/wssd_exon
{{
	/end1 exch def
	/start1 exch def
	/number1 exch def
	/y1 drop number1 mul def
	/y2 y1 5 add def
	/x1 start1 xstretch div cm def
	/x2 end1 xstretch div cm def
	newpath x1 y1  moveto x1 y2 lineto x2 y2 lineto x2 y1 lineto
	closepath
	fill stroke
}} def

/exon
{{
        /end1 exch def
        /start1 exch def
        /number1 exch def
        /y1 drop number1 mul def
        /y2 y1 10 add def
        /x1 start1 xstretch div cm def
        /x2 end1 xstretch div cm def
        newpath x1 y1  moveto x1 y2 lineto x2 y2 lineto x2 y1 lineto
        closepath
        fill stroke
}} def



% ------------------------------------------------
% DRAW MIROPEATS LINES

/drop {PAGE_LENGTH} {contig_count} div cm def
/xstretch {contig_stretch100} 100 div def

% draw the contig lines and names
{CTGS}

% draw the homology lines
/repcount {repeat_count} def
{SEGS}



% ------------------------------------------------
% DRAW REPEATMASKER
{RM}

% ------------------------------------------------
% DRAW DUPMASKER
{DM}

% ------------------------------------------------
% DRAW BED
{BED}

% ------------------------------------------------
% DRAW AXIS AND TICKS

{SCALE_WIDTH} cm setlinewidth 
0 {SCALE_BELOW} cm moveto
graphicmargin {SCALE_BELOW} cm lineto stroke

0 {SCALE_BELOW} cm moveto
0 {MIN_SCALE_STEP} {rightmost}{{
        /x1 exch {rightmost} div graphicmargin mul def
        newpath
        x1 {SCALE_BELOW} cm moveto
        x1 {SCALE_NUM_BELOW} cm lineto
        {TICK_WIDTH} cm setlinewidth
        stroke
        }} for

0 {SCALE_BELOW} cm moveto
0 {MAJ_SCALE_STEP} {rightmost}{{
        dup dup
        /x1 exch {rightmost} div graphicmargin mul def
        100 div truncate cvi cvr 10 div /numstring exch 10 string cvs def
        newpath
        x1 {SCALE_BELOW} cm moveto
        x1 {SCALE_NUM_BELOW} cm lineto
        currentpoint /y2 exch def /x2 exch def
        {TICK_WIDTH} cm setlinewidth
        stroke
        newpath
        x2 y2 moveto
        /{SCALE_FONT_TYPE} findfont {SCALE_FONT_SIZE} scalefont setfont
        numstring stringwidth pop -2 div -{SCALE_FONT_SIZE}  rmoveto
        numstring show

        }} for
        graphicmargin {SCALE_NUM_BELOW} cm -{SCALE_FONT_SIZE} add moveto (kbp) show

showpage
%%Trailer

"""


args.out.write(out)

