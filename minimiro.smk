import os 
import sys
import re

Q="" # query sequence
R="" # reference seqeunce 
REGIONS="" # region(s) to pull from reference (bed file)
SCORE=1000 # minimum peak dp socre 
SDIR=os.path.dirname(workflow.snakefile)
DEBUG=False

# delete if not debug
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))

rule all:
	input:
		pdf = "aln.pdf",

rule get_rgns:
	input:
		ref=R,
		fai=R + ".fai",
		rgn=REGIONS,
	output:
		fasta = tempd("temp/ref.fasta"),
		rgns = tempd("temp/one_based_rgns"),
	threads:1
	run:
		rgns = [ line.split()[0:3] for line in open(input["rgn"]).readlines() ]
		f = open(output["rgns"] , "w+")
		for contig, start, end in rgns:
			f.write( "{}:{}-{}\n".format( contig, int(start)+1, int(end) )  )
		f.close()
		shell("samtools faidx {input.ref} -r {output.rgns} > {output.fasta}")	
	

rule minimap2:
	input:
		ref = "temp/ref.fasta",
		query=Q,
	output:
		paf = tempd("temp/aln.paf"),
	threads: 32
	shell:"""
# YOU HAVE TO INCLUDE --cs FOR MINIMIRO TO WORK
minimap2 -x asm20 -s {SCORE} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}
"""


rule minimiro:
	input:
		paf = "temp/aln.paf",
	output:
		ps = tempd("temp/aln.ps"),
		pdf = tempd("temp/aln.pdf"),
	threads: 1
	shell:"""
{SDIR}/minimiro.py --paf {input.paf} \
	--bestn 1000 \
	-o {output.ps} && \
	ps2pdf {output.ps}
"""

rule move:
	input:
		pdf = "temp/aln.pdf",
	output:
		pdf = "aln.pdf",
	shell:"""
cp {input.pdf} {ouput.pdf}
"""




















