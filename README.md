# minimiro

## Install
```
git clone https://github.com/mrvollger/minimiro.git
```

## Usage: simple, no annotations
```
minimap2 -x asm20 -s {SCORE} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}

minimiro.py --paf {input.paf} --bestn 1000 -o {output.ps} && ps2pdf {output.ps}
```

## Usage: annotation pipeline 

## miropeats

### Configuration 
First make a file called `minimiro.yaml` that looks like this:
```
scores:  # minimum score threshold to plot, can specify multiple and it will make multiple pdfs
    - 50000
    - 150000


DEF_hg38_vs_CHM13: # run name for this comparison, will name output accordingly
    ref: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta
    regions:
        - chr8:6000000-13500000
    genes: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.refseq.genes.bed # bed12 gene file with ref coordiantes
    query: /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/t2t_rel3_glchr8/v8_hybrid/asm/HiCanu/chm13_hicanu_hifi_20k_glchr8v8.fa
    queryregions:
        - glchr8v8:6000000-13500000
    rc: False # if set to True the query will be reverse complemened before displaying
```

### Execution 
```
snakemake -j 200 -p -s /path/to/github/repo/minimiro.smk
```

If you don't want to clone the repo for you project you can use my local copy `/net/eichler/vol26/home/mvollger/projects/minimiro/minimiro.smk
` but use at your own risk because there will likely be modifications. 



### Example output
```
/net/eichler/vol26/projects/koren_hifi_asm/nobackups/SMN_hicanu/DEF_test/
```




## SEDEF
```
snakemake -k -p -s /path/to/repo/masker.smk -j 200 \
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree={resources.mem}G -pe serial {threads} -V -cwd -S /bin/bash -w n" \
    --drmaa-log-dir Masker/logs \
    sedef \
    --config fasta=$fasta sample=$sample
```

Here `$fasta` should point to an indexed fasta file, and `$sample` should be an informative sample name for the run.


