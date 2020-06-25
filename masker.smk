import os 
import sys
import re
import re
import pysam 
import pandas as pd
from datetime import date
today = date.today()
DATE =  today.strftime("%Y/%m/%d")

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
SDs="/net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.merged.mean.segdups.bed"
HG38="/net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta"
FOFN="/net/eichler/vol27/projects/sequence_data/nobackups/human/CHM13/PacBioHiFi/20_kbp_insert_hifi_beta/fastq.fofn"
SDADIR="/net/eichler/vol27/projects/hprc/nobackups/freeze1.0_assembly_eval/sda_collapse/software/centos7/SDA"

FASTA = os.path.abspath( config["fasta"] )
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"

# WILDCARDS
IDS = [ "{:08}".format(ID+1) for ID in range(len(open(FAI).readlines())) ]
SM = "asm"
if("sample" in config): SM = config["sample"]
SPECIES = "human"
if("species" in config): SPECIES = config["species"]
THREADS=8
if("threads" in config): THREADS = config["threads"]

wildcard_constraints:
	SM=SM,
	ID="\d+",

#
# OUTPUTS
#
HOSTNAME = os.getenv("HOSTNAME")
if("ocelot" in HOSTNAME and False):
	USER = os.getenv('USER')
	FASTA_FMT = f"/tmp/{USER}/rm_temp/{SM}_{{ID}}.fasta"
else:
	FASTA_FMT = f"temp/{SM}_{{ID}}.fasta"

DUP = f"{SM}_dupmasker.tbl"
COLOR = f"{SM}_dupmasker_colors.tbl"
RM = f"{SM}_repeatmasker.out"
RMBED = f"{SM}_repeatmasker.out.bed"
TRFBED = f"{SM}_trf.bed"
BED = f"{SM}_dupmasker_colors.bed"
DMHTML = f"{SM}_dupmasker_colors.html"
MASKED = f"{SM}_masked.fasta"
PAF = f"{SM}_to_hg38.paf"
COL = f"sda_out/coverage/{SM}.collapses.bed"
HTML = f"{SM}_sedef_out/SDs.browser.html"
DELTA = f"nucmer_out/{SM}.delta",

workdir: "Masker"

localrules: sda_collapse, run_sedef, sedef_browser, align_to_ref

#
# RULES
#
rule all:
	input:
		rm = RM,
		rmbed = RMBED,
		color = COLOR,
		dup = DUP,
		bed = BED, 
		dmhtml = DMHTML,
		masked = MASKED, 
		paf = PAF, 
		#collapse = COL,
		html = HTML, 
		#delta=DELTA, 

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
	threads: THREADS 
	shell:"""
echo "RM on {input.fasta}"
RepeatMasker \
	-s \
	-xsmall \
	-e ncbi \
	-species {SPECIES} \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta} 
	
if [ -f "{output.msk}" ]; then
    echo "maksed fasta exists"
else 
    echo "No repeats found, copying unmasked fasta to masked fasta"
	cp {input.fasta} {output.msk}
fi
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
	resources:
		mem=4,
	shell:"""
{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
"""

#
# trf rules 
#
rule trf:
	input:
		fasta = FASTA_FMT,
	output:
		dat = tempd(FASTA_FMT + ".dat")
	resources:
		mem=8,
	threads: 1
	shell:"""
trf {input.fasta} 2 7 7 80 10 50 15 -l 25 -h -ngs > {output.dat}
"""

rule trf_bed:
	input:
		dats = expand(rules.trf.output.dat, ID=IDS, SM=SM),
	output:
		bed = TRFBED,
	run:
		trf = []
		header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()
		for datf in input.dats:
			chrom = None
			sys.stderr.write( "\r" + datf )			
			with open(datf, 'r') as dat:
				for line in dat:
					splitline = line.split()
					if( line.startswith("Sequence:") ):
						chrom = int(line.split()[1].strip())
						#sys.stderr.write(chrom + "\n")
					elif( line.startswith("@") ):
						chrom = splitline[0][1:].strip() # grab everything after the @ in the first word
					else:
						# Catch index errors when line is blank
						try:
							# Check if in header sequence (all non-header lines start with an int: start pos)
							try:
								int(splitline[0])
							except ValueError:
								continue
							trf.append([chrom] + splitline[ 0: (len(header)-1) ] )
						except IndexError:
							pass
		trf = pd.DataFrame(trf, columns=header)
		print(trf.shape)
		
		trf["start"] = trf["start"].astype(int)
		trf.sort_values(by=["#chr", "start"], inplace=True)
		print("done sorting trf")

		trf.to_csv(output.bed, sep="\t", index=False)


#
# RepeatMasker merge output
#
rule mergeRM:
	input:
		outs = expand( rules.RepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		out=RM,
	resources:
		mem=4,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
"""

rule mergeRMbed:
	input:
		outs = expand( rules.RepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		bed=RMBED,
	resources:
		mem=4,
	threads: 1
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
# Make a masked fasta
#
rule rm_masked_fasta:
	input:
		msks = expand( rules.RepeatMasker.output.msk, ID=IDS, SM=SM),
	output:
		masked = temp(f"temp/{SM}.rm_masked.fasta"),
	resources:
		mem=8,
	threads:1
	shell:"""
cat {input.msks} > {output.masked}
"""

rule masked_fasta:
	input:
		fasta = rules.rm_masked_fasta.output.masked, 
		bed = TRFBED,
	output:
		fasta = MASKED
	resources:
		mem=8,
	threads:1
	shell:"""
seqtk seq -M {input.bed} {input.fasta} > {output.fasta}
"""
		


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
	resources:
		mem=4,
	threads: 1 
	run:
		f=open(output.dup, "w+")
		for dup in input.dups:
			f.write(open(dup).read())
		f.close()
	
		f=open(output.color, "w+")
		for idx, color in enumerate(input.colors):
			if(idx == 0):
				f.write(open(color).read())
			else:
				# all but the first line
				f.write( "".join( open(color).readlines()[1:] ) )
		f.close()
	


def hex_to_rgb(h):
	h = h.lstrip('#')
	return( ",".join( tuple( str(int(h[i:i+2], 16)) for i in (0, 2, 4))  )  )

rule DupMaskerBed:
	input:
		#color = rules.mergeDM.output.color,
		colors = expand(rules.DupMaskerColor.output.dupcolor,ID=IDS, SM=SM),
	output:
		bed =  BED
	resources:
		mem=4,
	threads: 1 
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

		color.sort_values(by=["chr", "chrStart"], inplace=True)
		out = color[ ["chr", "chrStart", "chrEnd", "Repeat", "score", "strand", "chrStart", "chrEnd", "rgb"] ]
		out.to_csv(output["bed"], sep="\t", header=False, index=False)			

rule DupMaskerHTML:
	input:
		bed=BED,
	output:
		html = DMHTML,
	resources:
		mem=4,
	threads: 1 
	run:
		html = open(f"{SDIR}/templates/dupmasker.html").read()
		open(output["html"], "w+").write(html.format(DATE=DATE, SM=SM))

			
		
#
# nucmer rules
#
rule nucmer:    
	input:
		q = FASTA,
		r = HG38,
	output:
		delta = DELTA
	threads: 128
	resources:
		mem=4,
	shell:"""
mkdir -p nucmer_out
module load mummer/4.0beta
nucmer -t {threads} -l 100 --maxmatch  --prefix=nucmer_out/{SM} {input.r} {input.q}
"""



rule tiling:
	input:
		delta = rules.nucmer.output.delta,
	output:
		tiling = str(rules.nucmer.output.delta)+ ".tiling"		
	threads: 1
	resources:
		mem=8,
	shell:"""
module load mummer/4.0beta
show-tiling {input.delta} > {output.tiling}
"""

rule diffs:	
	input:
		delta = rules.nucmer.output.delta,
	output:
		rdiff = str(rules.nucmer.output.delta)+ ".r.diff",
		qdiff = str(rules.nucmer.output.delta)+ ".q.diff",
	threads: 1
	resources:
		mem=8,
	shell:"""
module load mummer/4.0beta
show-diff -H -r {input.delta} > {output.rdiff}
show-diff -H -q {input.delta} > {output.qdiff}
"""




rule dnadiff:	
	input:
		delta = rules.nucmer.output.delta,
	output:
	threads: 1
	resources:
		mem=8,
	shell:"""
module load mummer/4.0beta
dnadiff -p nucmer_out/{SM} -d {input.delta} 
"""


#
# sedef rules
#
rule run_sedef:
	input:
		fasta = MASKED,
	output:
		tmpf = temp(f"/tmp/mvollger/sedef/sedef_{SM}.fasta"),
		bed = f"{SM}_sedef_out/final.bed"
	resources:
		mem=8,
	threads: 80
	shell:"""
module load boost/1.68.0
PATH=/net/eichler/vol26/projects/koren_hifi_asm/nobackups/software/sedef:$PATH
rm -f {output.tmpf} {output.tmpf}.fai
sedef.sh -f -o $( dirname {output.bed} ) -j {threads} -t {output.tmpf} {input.fasta}
"""


html=f"""<h2>Description</h2>
A track showing segmental duplications as generated by <a href="https://github.com/vpc-ccg/sedef">sedef</a> on the assembly: {SM}.
Column descriptions can be found here: <a href="https://github.com/vpc-ccg/sedef#output">columns</a>.
<p></p>

<h2>Methods</h2>
<p>
Sedef was run with default parameters on a RepeatMasked genome assembly.
The resulting output (final.bed) was then converted into a browser friendly format (bed9 + extra fields).
</p>

<h2>Display Conventions and Configuration</h2>
<p>
Purple: less than 90% similarity
<br>Light to dark gray: 90 - 98% similarity
<br>Yellow: 98 - 99% similarity
<br>Orange: greater than 99% similarity
</p>


<h2>Data access</h2>
<!-- Optional, links for downloading the data if there is a separate source beside the
  browser  -->
<p></p>

<h2>Release history</h2>
<ol>
  <li> {SM}
</ol>

<h2>Contact</h2>
<p>Contact <a href="mailto:mvollger@uw.edu"> Mitchell R. Vollger &lt;mvollger@uw.edu&gt;</a>
</p>

<h2>Credits</h2>
<p></p>
"""
	
rule sedef_browser:
	input:
		bed = rules.run_sedef.output.bed,
	output:
		bed = f"{SM}_sedef_out/SDs.browser.bed",
		lowid = f"{SM}_sedef_out/SDs.lowid.browser.bed",
		#html = f"{SM}_sedef_out/SDs.browser.html",
		html = HTML,
	resources:
		mem=8,
	threads: 1
	run:
		html = open(f"{SDIR}/templates/sedef.html").read()
		open(output["html"], "w+").write(html.format(DATE=DATE, SM=SM))
		shell("{SDIR}/scripts/sedef_to_bed.py {input.bed} {output.bed} {output.lowid}")

rule sedef:
	input:
		rules.sedef_browser.output,





#
# collapse anlysis 
#	
rule sda_collapse:
	input:
		fasta = FASTA,
		fofn = FOFN,
	output:
		collapse = COL,
	resources:
		mem=8,
	threads: 1
	shell:"""
{SDADIR}/SDA denovo \
        --drmaa " -l mfree={{resources.mem}}G -pe serial {{threads}} -l h_rt=128:00:00 -V -cwd -S /bin/bash " \
        --threads 400 \
        --platform ccs \
        --fofn {input.fofn} --ref {input.fasta} --pre {SM}  \
        {output.collapse}
"""



#
# minimap2 rules
#
rule align_to_ref:
	input:
		q = FASTA,
		r = HG38,
	output:
		paf = PAF,
	threads: 128
	resources:
		mem=4,
	shell:"""
minimap2 -I 8G -2K 2000m -t {threads} \
	--secondary=no -c --eqx \
	-x asm20 -s 200000 \
	 -z 10000 -r 50000 \
	--paf-no-hit \
	{input.r}  {input.q} > {output.paf}
"""

rule split_align:
	input:
		q = FASTA,
		r = HG38,
	output:
		split = temp(f"split_alignments/{SM}.split.fasta"),
		bam = f"split_alignments/{SM}.split.bam",
		bai = f"split_alignments/{SM}.split.bam.bai",
	threads: 128
	resources:
		mem=4,
	shell:"""
~/projects/utility/split_fasta.py -n 5000 {input.q} > {output.split}

minimap2 -I 8G -2K 2000m -t {threads} \
	--secondary=no -c --eqx -Y \
	-ax asm20  \
	-r 50000 \
	{input.r}  {output.split} | samtools view -b - | samtools sort - > {output.bam}

samtools index {output.bam}

"""

rule split_bed:
	input:
		bam = rules.split_align.output.bam,
		bai = rules.split_align.output.bai,
	output:
		bed = f"split_alignments/{SM}.split.bed",
	threads: 1
	resources:
		mem=8,
	shell:"""
~/projects/hifi_asm/scripts/samIdentity.py --bed {input.bam} > {output.bed}
"""

FLANK=50000
MERGE_DIST = max(50000, FLANK+1)
MINSIZE=MERGE_DIST
rule inter_sd:
	input:
		sd=SDs,
		bed = rules.split_bed.output.bed,
		fai = HG38 + ".fai",
	output:
		bed = f"split_alignments/{SM}.split.sd.bed",
		nosd = f"split_alignments/{SM}.split.nosd.bed",
		flank = f"split_alignments/{SM}.flank.bed",
		sdflank = f"split_alignments/{SM}.split.sdflank.bed",
	threads: 1
	resources:
		mem=8,
	shell:"""
bedtools intersect -header -f 0.5 -wa -a {input.bed} -b {input.sd} > {output.bed}
bedtools intersect -header -v  -a {input.bed} -b {input.sd} > {output.nosd}

bedtools sort -i {input.sd} | bedtools merge -d {MERGE_DIST} -i - | \
	awk '$3-$2 > {MINSIZE}{{print $0}}' | \
	bedtools flank -g {input.fai} -b {FLANK} -i - > {output.flank}

bedtools intersect -header -f 0.5 -wa -a {input.bed} -b {output.flank} > {output.sdflank}

"""

rule sd_genes:
	input:	
		bed = rules.inter_sd.output.bed,
	output:
		bed = f"split_alignments/{SM}.split.sd.genes.bed",
	resources:
		mem=8,
	threads: 1
	shell:"""
bedtools intersect -wa -wb -a ~/assemblies/hg38/hg38.gene.locations.merged.bed -b {input.bed} | bedtools groupby -g 1,2,3,4 -c 12 -o mean > {output.bed}
"""

rule hg38:
	input:
		paf=PAF,
		bam = rules.split_align.output.bam,
		bai = rules.split_align.output.bai,
		bed = rules.split_bed.output.bed,
		sd = rules.inter_sd.output,
		genes = rules.sd_genes.output,



rule new_sd:
	input:
		sedef = rules.sedef_browser.output.bed,
		fasta = FASTA,
		fai = FASTA + ".fai",
	output:
		sdw = f"split_alignments/{SM}.sd.windows.bed",
		fasta = f"split_alignments/{SM}.sd.windows.fasta",
	threads: 1
	resources:
		mem=8,
	shell:"""
bedtools makewindows -g {input.fai} -w 5000  | bedtools intersect -wb -a {input.sedef} -b - > {output.sdw}
bedtools getfasta -fi {input.fasta} -bed {output.sdw} > {output.fasta}
"""
