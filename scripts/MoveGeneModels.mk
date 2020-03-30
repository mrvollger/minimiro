REF=/net/eichler/vol26/home/mvollger/projects/minimiro/data/hg38.fa
GENES=/net/eichler/vol26/home/mvollger/projects/minimiro/data/hg38.genes.bed  # must be bed12 format with genes 
QUERY=/net/eichler/vol26/projects/koren_hifi_asm/nobackups/koren_asm_2020_01_06/chm13_hicanu_hifi_20k.fasta
THREADS=64
MAXGB=64

all: \
	query.genes.12.bed \
	query.genes.6.bed

ref.cds.fasta:
	bedtools getfasta -name -split -fi $(REF) -bed $(GENES) > $@
	
aln.bam: ref.cds.fasta
	minimap2 --eqx -ax splice -C5 -O6,24 -B4 -t $(THREADS) --cap-sw-mem=$(MAXGB)g \
		$(QUERY) $< | \
		samtools view -b - | \
		samtools sort -m $(MAXGB)G - > $@ && samtools index $@

query.genes.12.bed: aln.bam
	bedtools bamtobed -bed12 -i $< > $@

query.genes.6.bed: query.genes.12.bed
	bedtools bed12tobed6 -i $< > $@

