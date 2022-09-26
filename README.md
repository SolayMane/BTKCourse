# BTKCourse
## BlobToolKit (v3.3.0) instaltion 
pip install blobtoolkit


This reposotory contains all codes used to analyse the Botrytis genome assemblies
## The rawd reads were trimmed suing sickle
## The trimmed reads were then assembled using clc genomics
## The 112-name assembly was assesed using blobtoolkit --> suspected contaminatioinn from bacetria!!
## I have used kraken to extracte only the Ascomycota reads and taxid children 
````bash
kraken2 --db /sanhome2/Comparative_genomics/Vcovid19/kraken-VBE-DB/ \
--threads 56 \
--report kraken/ngs.report  \
--paired trimmed/112-name.trim.R1.fastq.gz trimmed/112-name.trim.R2.fastq.gz >> kraken/112-name.kraken

extract_kraken_reads.py -k 112-name.kraken --include-children -t 4890 \
-s ../Assembly/trimmed/112-name.trim.R1.fastq.gz \
-s2 ../Assembly/trimmed/112-name.trim.R2.fastq.gz \
-o 112-name.trim.asco.R1.fastq.gz \
-o2 112-name.trim.asco.R2.fastq.gz \
-r ngs.report \
--max 10000000000 \
--fastq-output

````

## I reassembled the extracted reads to asses the assembly using blobtoolkit

## Adding data to a dataset (Create blobtoolkit project)
### Create a metadata file
1. Create yaml file to describe the project:
````bash
vim ~/fistrAssembly.yaml
````
````bash
assembly:
  accession: GCA_001028725.1
  alias: S_venezuelensis_HH1
  bioproject: PRJEB530
  biosample: AMD00012916
  record_type: contig
taxon:
  name: Strongyloides venezuelensis
````
2. Run the following command:
````bash
blobtools create --fasta /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa \
--meta asm.yaml  --taxid 75913 \
--taxdump /home1/software/blobtoolkit/taxdump Deconta_asm/
````


### Adding hits
1. Run the blastn: 
````bash
blastn -db /home1/software/blobtoolkit/nt/nt \
-query /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa \
-outfmt "6 qseqid staxids bitscore std" \
-max_target_seqs 10 -max_hsps 1 \
-evalue 1e-25 -num_threads 56 \
-out asm.ncbi.blastn.out
````

2. Run diamond:
````bash
diamond blastx --query /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa\
--db /home1/software/blobtoolkit/uniprot/reference_proteomes.dmnd \
--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
--sensitive --max-target-seqs 1 \
--evalue 1e-25 \
--threads 56 > asm.diamond.blastx.out
````

3. Add hits with blobtools command:
````bash
blobtools add --hits asm.ncbi.blastn.out \
--hits asm.diamond.blastx.out \
--taxrule bestsumorder \
--taxdump /home1/software/blobtoolkit/taxdump Deconta_asm/

### Adding Coverage
1. Run minimap:
````bash
minimap2 -ax sr -t 56 /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa \
../kraken/112-name.trim.asco.R1.fastq ../kraken/112-name.trim.asco.R2.fastq \
| samtools sort -@56 -O BAM -o asm.DRR008460.bam -
````
2. add coverage using blobtools command:
````bash
blobtools add --cov asm.bam --threads 56 Deconta_asm/
```` 
### Adding buscos
1. Run busco on the genome assemlby 
````bash
busco -i /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa \
-l helotiales_odb10 -o 112-name.asco -m genome --cpu 56
````
2. Add busco file using blobtolls command :
````bash
blobtools add --busco busco/112-name.asco/run_helotiales_odb10/full_table.tsv Deconta_asm/
````
### View the project on blobtool viewer
1. Run the following command to initialize the viewer
````bash
blobtools view --remote Deconta_asm/
````
2. Open your browser and go the URL indicated in the previous command (e.g http://localhost:8001/view/Deconta_asm)
