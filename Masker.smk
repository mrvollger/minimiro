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

wildcard_constraints:
	SM=SM,
	ID="\d+",

#
# OUTPUTS
#
FASTA_FMT = f"temp/{SM}_{{ID}}.fasta"
DUP = f"{SM}_dupmasker.tbl"
COLOR = f"{SM}_dupmasker_colors.tbl"
RM = f"{SM}_repeatmasker.out"
BED = f"{SM}_dupmasker_colors.bed"

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
		#cat = tempd(FASTA_FMT + ".cat.gz"),
		#tbl = tempd(FASTA_FMT + ".tbl"),
		#msk = tempd(FASTA_FMT + ".masked"),
	resources:
		mem=8,
	threads:8
	shell:"""
RepeatMasker \
	-e ncbi \
	-species {SPECIES} \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta}
"""


rule DupMasker:
	input:
		fasta = FASTA_FMT,
		out = rules.RepeatMasker.output.out,
	output:
		dups = tempd(FASTA_FMT + ".duplicons"),
		#dupout = tempd(FASTA_FMT + ".dupout"),
	resources:
		mem=8,
	threads:8
	shell:"""
{SDIR}/bin/DupMaskerParallel \
	-engine ncbi \
	{input.fasta}
"""
"""
DupMasker \
"""


rule DupMaskerColor:
	input:
		dups = rules.DupMasker.output.dups,
	output:
		dupcolor = tempd(FASTA_FMT + ".duplicons.extra"),
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""


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


rule mergeRM:
	input:
		outs = expand( rules.RepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		out=RM,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
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
			color = pd.read_csv(input.color, sep="\s+")
			colors.append(color)

		color = pd.concat(colors, ignore_index=True)
		color.sort_values(by=["chr", "chrStart"], inplace=True)

		color["strand"] = "+"
		color["strand"][ color["orient"] == "R" ] = "-"

		color["rgb"] = color["color"].map(hex_to_rgb)
		color["score"] = 0

		out = color[ ["chr", "chrStart", "chrEnd", "Repeat", "score", "strand", "rgb"] ]

		print(color)
			
			

			
		

		




