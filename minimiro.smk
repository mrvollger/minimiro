import os 
import sys
import re
import re
import pandas as pd 

configfile: "minimiro.yaml"
SDIR=os.path.dirname(workflow.snakefile)
DEBUG=True


SMS = list(config.keys())
SEQS=["ref", "query"]
SCORES=config["scores"];  SMS.remove("scores") # minimum peak dp socre 

RS = {}
RGNS = {}
QS = {}
QRGNS = {}
RCS = {}
GENES = {}
for SM in SMS:
	RS[SM] = config[SM]["ref"] # reference seqeunce 
	assert os.path.exists(RS[SM]+".fai")
	RGNS[SM]	=	config[SM]["regions"] # region(s) to pull from reference 
	QS[SM]		=	config[SM]["query"] # query sequence
	assert os.path.exists(QS[SM]+".fai")
	QRGNS[SM]	=	config[SM]["queryregions"] # region(s) to pull from query
	if("rc" in config[SM]):
		RCS[SM] = config[SM]["rc"] 
	else: 
		RCS[SM] = False
	GENES[SM] = config[SM]["genes"]	

wildcard_constraints:
	SEQ="|".join(SEQS),
	SM="|".join(SMS),
	SCORE="\d+",

rule all:
	input:
		pdf	= expand("minimiro_smk_out/{SM}_{SCORE}_aln.pdf", SM=SMS, SCORE=SCORES),

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
module load RepeatMasker/4.1.0
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
	threads:1
	shell:"""
module load RepeatMasker/4.1.0
DupMasker -engine ncbi \
	{input.fasta}
"""
#-pa {threads} \

rule DupMaskerColor:
	input:
		dups = rules.DupMasker.output.dups,
	output:
		dupcolor = "temp/{SM}_{SEQ}.fasta.duplicons.extra",
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""

#
# rules to get genes
#
def get_genes(wildcards):
	SM = str(wildcards.SM)
	return(GENES[SM])

def get_ref_bed(wildcards):
	SM = str(wildcards.SM)
	rtn = []
	for rgn in RGNS[SM]:
		match = re.match("(.+):(\d+)-(\d+)", rgn.strip())
		if(match):
			rtn.append('{}\t{}\t{}\n'.format(*match.groups()))
		else:
			rtn.append('{}\t{}\t{}\n'.format(rgn, 0, 1000000000))
	return( rtn )

rule get_cds:
	input:
		fasta = get_ref,
		bed = get_genes,
	output:
		#fasta = "temp/{SM}.genes.fasta",
		bed = "temp/{SM}.ref.genes.bed",
		tmp = temp("temp/tmp.{SM}.ref.genes.bed"),
	params:
		bed = get_ref_bed,
	run:
		# read in regular symbol names 
		convert = {}
		for line in open(f"{SDIR}/data/gene_conversion.txt").readlines():
			t = line.strip().split()
			if(len(t) > 1):
				convert[  t[1].split(".")[0] ] = t[0]
	
		# make a bed file with all the corrdiantes in the correct space
		rtn = ""
		for bed in params["bed"]:
			shell('bedtools intersect -a {input.bed} -b <(printf "{bed}") | bedtools bed12tobed6 -i /dev/stdin > {output.tmp}')
			chrm, start, end = bed.split(); start = int(start);
			name = f"{chrm}:{start}-{end}"
			for line in open(output["tmp"]).readlines():
				t = line.strip().split()
				t[0] = name
				t[1] = int(t[1]) - start	
				t[2] = int(t[2]) - start	
				gene_id = t[3].split(".")[0]
				if(gene_id in convert): 
					t[3] = convert[gene_id]
					rtn += "{}\t{}\t{}\t{}\t{}\t{}\n".format(*t)
		open(output["bed"], "w+").write(rtn)
		
		# get all the cds seqeunces to map to the query 
		#shell("bedtools getfasta -name -split -fi {input.fasta} -bed {input.bed} > {output.fasta}")

#
# make the alignments and the miropeats 	
#
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
minimap2 -x asm20 -r 200000 -s {params.score} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}
"""


rule minimiro:
	input:
		paf = rules.minimap2.output.paf,
		rmout = expand("temp/{{SM}}_{SEQ}.fasta.out", SEQ=SEQS),
		dmout = expand("temp/{{SM}}_{SEQ}.fasta.duplicons.extra", SEQ=SEQS),
		genes = rules.get_cds.output.bed,
	output:
		ps	= "minimiro_smk_out/{SM}_{SCORE}_aln.ps",
		pdf	= "minimiro_smk_out/{SM}_{SCORE}_aln.pdf",
	threads: 1
	shell:"""
{SDIR}/minimiro.py --paf {input.paf} \
	--rm {input.rmout} \
	--dm {input.dmout} \
	--bed {input.genes} \
	--bestn 1000 \
	-o {output.ps} && \
	ps2pdf {output.ps} {output.pdf}
"""










