import os 
import sys
import re

configfile: "minimiro.yaml"
SDIR=os.path.dirname(workflow.snakefile)
DEBUG=True


SMS = list(config.keys())

SCORES=config["scores"];  SMS.remove("scores") # minimum peak dp socre 

RS = {}
RGNS = {}
QS = {}
QRGNS = {}
RCS = {}
for SM in SMS:
	RS[SM]		=	config[SM]["ref"] # reference seqeunce 
	assert os.path.exists(RS[SM]+".fai")
	
	RGNS[SM]	=	config[SM]["regions"] # region(s) to pull from reference 
	
	QS[SM]		=	config[SM]["query"] # query sequence
	assert os.path.exists(QS[SM]+".fai")
	
	QRGNS[SM]	=	config[SM]["queryregions"] # region(s) to pull from query
	if("rc" in config[SM]):
		RCS[SM] = config[SM]["rc"] 
	else: 
		RCS[SM] = False


SEQS=["ref", "query"]
wildcard_constraints:
	SEQ="|".join(SEQS),
	SM="|".join(SMS),
	SCORE="\d+",


# delete if not debug
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))

def get_ref(wildcards):
	SM = str(wildcards.SM)
	return(RS[SM])
def get_query(wildcards):
	SM = str(wildcards.SM)
	return(QS[SM])

def get_ref_rgns(wildcards):
	SM = str(wildcards.SM)
	return( " ".join( RGNS[SM] ) )
def get_query_rgns(wildcards):
	SM = str(wildcards.SM)
	return( " ".join( QRGNS[SM] ) )

def get_rc(wildcards):
	SM = str(wildcards.SM)
	return(RCS[SM])

rule get_rgns:
	input:
		ref=get_ref, 
		query=get_query,
	output:
		ref = tempd("temp/{SM}_ref.fasta"),
		query = tempd("temp/{SM}_query.fasta"),
	params:
		rgns = get_ref_rgns,
		qrgns = get_query_rgns,
		rc = get_rc,
	threads:1
	run:
		shell("samtools faidx {input.ref} {params.rgns} > {output.ref}")	
		if(params["rc"]):
			shell("samtools faidx {input.query} {params.qrgns} | seqtk seq -r - > {output.query}")	
		else:	
			shell("samtools faidx {input.query} {params.qrgns}- > {output.query}")	

rule RepeatMasker:
	input:
		fasta = "temp/{SM}_{SEQ}.fasta",
	output:
		out = tempd("temp/{SM}_{SEQ}.fasta.out"),
		cat = tempd("temp/{SM}_{SEQ}.fasta.cat"),
		tbl = tempd("temp/{SM}_{SEQ}.fasta.tbl"),
		msk = tempd("temp/{SM}_{SEQ}.fasta.masked"),
	resources:
		mem=8,
	threads:8
	shell:"""
module load perl/5.14.2 RepeatMasker/3.3.0
RepeatMasker \
	-species human \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta}
"""


rule DupMasker:
	input:
		fasta = "temp/{SM}_{SEQ}.fasta",
		out = rules.RepeatMasker.output.out,
	output:
		dups = "temp/{SM}_{SEQ}.fasta.duplicons",
	shell:"""
module load perl/5.14.2 RepeatMasker/3.3.0
DupMasker {input.fasta}
"""

rule DupMaskerColor:
	input:
		dups = rules.DupMasker.output.dups,
	output:
		dupcolor = "temp/{SM}_{SEQ}.fasta.duplicons.extra",
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""
		
 
def get_score(wildcards):
	return( int(str(wildcards.SCORE)))

rule minimap2:
	input:
		ref = rules.get_rgns.output.ref,
		query = rules.get_rgns.output.query,
	output:
		paf = tempd("temp/{SM}_{SCORE}_aln.paf"),
	params:
		score = get_score,
	threads: 16
	shell:"""
# YOU HAVE TO INCLUDE --cs FOR MINIMIRO TO WORK
minimap2 -x asm20 -s {params.score} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}
"""


rule minimiro:
	input:
		paf = rules.minimap2.output.paf,
		rmout = expand("temp/{{SM}}_{SEQ}.fasta.out", SEQ=SEQS),
		dmout = expand("temp/{{SM}}_{SEQ}.fasta.duplicons.extra", SEQ=SEQS),
	output:
		ps	= "minimiro_smk_out/{SM}_{SCORE}_aln.ps",
		pdf	= "minimiro_smk_out/{SM}_{SCORE}_aln.pdf",
	threads: 1
	shell:"""
{SDIR}/minimiro.py --paf {input.paf} \
	--rm {input.rmout} \
	--dm {input.dmout} \
	--bestn 1000 \
	-o {output.ps} && \
	ps2pdf {output.ps} {output.pdf}
"""



rule all:
	input:
		pdf	= expand(rules.minimiro.output.pdf, SM=SMS, SCORE=SCORES),



















