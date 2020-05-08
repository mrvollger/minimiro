import os 
import sys
import re
import re
import pysam 
import pandas as pd

SDIR=os.path.dirname(workflow.snakefile)
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

# delete if not debug
DEBUG=True
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))

#
# INPUTS
#
FASTA = os.path.abspath( config["fasta"] )
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"

# WILDCARDS
IDS = [ "{:08}".format(ID+1) for ID in range(len(open(FAI).readlines())) ]
SM = "asm"
if("sample" in config): SM = config["sample"]
SPECIES = "human"
if("species" in config): SPECIES = config["species"]
THREADS=16
if("threads" in config): THREADS = config["threads"]

wildcard_constraints:
	SM=SM,
	ID="\d+",

#
# OUTPUTS
#
HOSTNAME = os.getenv("HOSTNAME")
if("ocelot" in HOSTNAME):
	USER = os.getenv('USER')
	FASTA_FMT = f"/tmp/{USER}/rm_temp/{SM}_{{ID}}.fasta"
else:
	FASTA_FMT = f"temp/{SM}_{{ID}}.fasta"

DUP = f"{SM}_dupmasker.tbl"
COLOR = f"{SM}_dupmasker_colors.tbl"
RM = f"{SM}_repeatmasker.out"
RMBED = f"{SM}_repeatmasker.out.bed"
BED = f"{SM}_dupmasker_colors.bed"
MASKED = f"{SM}_rm_masked.fasta"

workdir: "Masker"


#
# RULES
#
rule all:
	input:
		rm = RM,
		color = COLOR,
		dup = DUP,
		bed = BED, 
		masked = MASKED, 

rule split_fasta:
	input:
		fasta = FASTA,
	output:
		fastas = tempd(expand(FASTA_FMT, ID=IDS)),
	threads: 1
	resources:
		mem=8
	run:
		fasta = pysam.FastaFile(input["fasta"])
		for out, name in zip(output["fastas"], fasta.references):
			seq = fasta.fetch(name)
			open(out, "w+").write( ">{}\n{}\n".format(name, seq) )


rule RepeatMasker:
	input:
		fasta = FASTA_FMT,
	output:
		out = tempd(FASTA_FMT + ".out"),
		msk = tempd(FASTA_FMT + ".masked"),
		#cat = tempd(FASTA_FMT + ".cat.gz"),
		#tbl = tempd(FASTA_FMT + ".tbl"),
	resources:
		mem=8,
	log:
		f"logs/RM.{SM}.{{ID}}.log",
	threads: int( THREADS / 2 )
	shell:"""
echo "RM on {input.fasta}"
RepeatMasker \
	-xsmall \
	-e ncbi \
	-species {SPECIES} \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta} 
	

#-frag 10000000 \
"""



rule DupMaskerRM:
	input:
		fasta = FASTA_FMT,
		out = rules.RepeatMasker.output.out,
	output:
		dupout = tempd(FASTA_FMT + ".dupout"),
	resources:
		mem=8,
	threads: THREADS
	shell:"""
{SDIR}/bin/DupMaskerParallel \
	-pa {threads} -dupout \
	-engine ncbi \
	{input.fasta}
"""

rule DupMasker:
	input:
		fasta = FASTA_FMT,
		out = rules.RepeatMasker.output.out,
		dupout = rules.DupMaskerRM.output.dupout,
	output:
		dups = tempd(FASTA_FMT + ".duplicons"),
	resources:
		mem=8,
	threads:1
	shell:"""
{SDIR}/bin/DupMaskerParallel \
	-engine ncbi \
	{input.fasta}
"""


rule DupMaskerColor:
	input:
		dups = rules.DupMasker.output.dups,
	output:
		dupcolor = tempd(FASTA_FMT + ".duplicons.extra"),
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""





#
# RepeatMasker merge output
#
rule mergeRM:
	input:
		outs = expand( rules.RepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		out=RM,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
"""

rule masked_fasta:
	input:
		msks = expand( rules.RepeatMasker.output.msk, ID=IDS, SM=SM),
	output:
		masked = MASKED,
	threads:1
	shell:"""
cat {input.msks} > {output.masked}
"""


rule mergeRMbed:
	input:
		outs = expand( rules.RepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		bed=RMBED,
	run:
		rms = []
		for frm in input.outs:
			sys.stderr.write( "\r" + frm )
			rms.append( pd.read_csv(frm, delim_whitespace=True, header=None, skiprows=[0,1,2], comment="*",
				names = ["score", "div", "del", "ins", "#contig", "start", "end",
					 "q_left", "strand", "repeat", "class", "r_st", "r_en", "r_left", "id"]) )
			
		rm = pd.concat(rms, ignore_index=True)
		
		rm[['repeat_class','repeat_subclass']] = rm["class"].str.split("/", expand=True)
		#print(rm[['repeat_class','repeat_subclass']])	
		rm.repeat_subclass.fillna(value=".", inplace=True)
		
		rm.loc[rm.strand=="C", "strand"] = "-"

		outcolumns = ["#contig", "start", "end", "repeat", "score", "strand", "repeat_class", 
			"repeat_subclass", "div", "del", "ins", "q_left", "r_st", "r_en", "r_left"]
		
		rm = rm[outcolumns]	
		rm.sort_values(by=["#contig", "start"], inplace=True)
		rm.to_csv(output.bed, sep="\t", index=False)

#
# DupMakser merge output
#

rule mergeDM:
	input:
		dups = expand(rules.DupMasker.output.dups,ID=IDS, SM=SM),
		colors = expand(rules.DupMaskerColor.output.dupcolor,ID=IDS, SM=SM),
	output:
		dup = DUP,
		color = COLOR,
	shell:"""
cat {input.dups} > {output.dup}

head -n 1 {input.colors[0]} > {output.color} && \
	tail -q -n +2 {input.colors} >> {output.color}
"""


def hex_to_rgb(h):
	h = h.lstrip('#')
	return( ",".join( tuple( str(int(h[i:i+2], 16)) for i in (0, 2, 4))  )  )

rule DupMaskerBed:
	input:
		#color = rules.mergeDM.output.color,
		colors = expand(rules.DupMaskerColor.output.dupcolor,ID=IDS, SM=SM),
	output:
		bed =  BED
	run:
		colors = []
		for f in input.colors:
			#chr    chrStart	chrEnd  orient  Repeat  color   width   offset
			color = pd.read_csv(f, sep="\s+")
			colors.append(color)

		color = pd.concat(colors, ignore_index=True)

		# these rows are no good, they come from contigs that have messed up results. fix TODO
		bad= (color["chr"] == 0 ) & ( color["chrEnd"] == "#BEBEBE")
		color.drop(color[bad].index, inplace = True)

		color.sort_values(by=["chr", "chrStart"], inplace=True)

		color["strand"] = "+"
		color.loc[color["orient"] == "R", "strand"] = "-"

		color["rgb"] = color["color"].map(hex_to_rgb)
		color["score"] = 0

		out = color[ ["chr", "chrStart", "chrEnd", "Repeat", "score", "strand", "rgb"] ]

		out.to_csv(output["bed"], sep="\t", header=False, index=False)			

			
		

		




