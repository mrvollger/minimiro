import os 
import sys
import re
import re
import pysam 
SDIR=os.path.dirname(workflow.snakefile)
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

# delete if not debug
DEBUG=True
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))


FASTA = os.path.abspath( config["fasta"] )
FAI = FASTA + ".fai"
SPECIES = "human"
ID_FMT = "{:08}"
IDS = [ ID_FMT.format(ID+1) for ID in range(len(open(FAI).readlines())) ]

workdir: "Masker"

# OUTPUTS
FASTA_FMT = "temp/{ID}.fasta"
DUP = "dupmasker.tbl"
COLOR = "dupmasker_colors.tbl"
RM = "repeatmasker.out"

assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"



rule all:
	input:
		rm = RM,
		color = COLOR,
		dup = DUP,

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
		cat = tempd(FASTA_FMT + ".cat.gz"),
		tbl = tempd(FASTA_FMT + ".tbl"),
		msk = tempd(FASTA_FMT + ".masked"),
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
		dups = FASTA_FMT + ".duplicons",
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
		dupcolor = FASTA_FMT + ".duplicons.extra",
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""


rule mergeDM:
	input:
		dups = expand(rules.DupMasker.output.dups,ID=IDS),
		colors = expand(rules.DupMaskerColor.output.dupcolor,ID=IDS),
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
		outs = expand( rules.RepeatMasker.output.out, ID=IDS),
	output:
		out=RM,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
"""




