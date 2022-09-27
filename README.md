# BTKCourse
## BlobToolKit (v3.3.0) instaltion 
<code>pip install blobtoolkit</code>

## Demo on BTK
```mermaid
graph TB
subgraph CASE 2
    er(Raw Reads) -->|extracte target reads| aa(Newdenovo assembly)
    aa -- <i>Blobtoolkit --> t(Clean Assembly)
   
  end
  subgraph CASE 1
    tt(Raw Reads) --> a
    a[denovo assembly] -- <i>Blobtoolkit --> b(Contaminated assembly!)

    a --> zz(Clean Assembly!)
    b -->|filter| EE(target group contigs)
  end
  
        classDef green fill:#93FF33,stroke:#333,stroke-width:2px
        classDef blue fill:#00FA9A,stroke:#333,stroke-width:4px
       
        class g,a,h green
        class b,c,d,e,f blue
 ```      

## Adding data to a dataset (Create blobtoolkit project)
### 1. Create a metadata file
1. Create a project directory
````bash
sudo mkdir ~/btk/
````
2. Create yaml file to describe the project:
````bash
vim ~/fistrAssembly.yaml
````
````bash
assembly:
  alias: B_cinera_112
  record_type: contig
taxon:
  name: Botrytis cinerea
````
3. Run the following command to creat btk directory project
````bash
blobtools create --fasta ~/mygenome.fasta \
--meta ~/fistrAssembly.yaml  --taxid 75913 \
--taxdump /home1/software/blobtoolkit/taxdump ~/btk
````

### 2. Adding hits
1. Run the blastn: 
````bash
blastn -db /home1/software/blobtoolkit/nt/nt \
-query ~/mygenome.fasta \
-outfmt "6 qseqid staxids bitscore std" \
-max_target_seqs 10 -max_hsps 1 \
-evalue 1e-25 -num_threads 56 \
-out asm.ncbi.blastn.out
````

2. Run diamond:
````bash
diamond blastx --query ~/mygenome.fasta \
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
````
### 3. Adding Coverage
1. Run minimap:
````bash
minimap2 -ax sr -t 56 /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa \
../kraken/112-name.trim.asco.R1.fastq ../kraken/112-name.trim.asco.R2.fastq \
| samtools sort -@56 -O BAM -o asm.DRR008460.bam -
````
2. add coverage using blobtools command:
````bash
blobtools add --cov asm.bam --threads 56 ~/btk
```` 
### 4. Adding buscos
1. Run busco on the genome assemlby 
````bash
busco -i ~/mygenome.fasta \
-l helotiales_odb10 -o 112-name.asco -m genome --cpu 56
````
2. Add busco file using blobtolls command :
````bash
blobtools add --busco busco/112-name.asco/run_helotiales_odb10/full_table.tsv ~/btk
````
### 5. View the project on blobtool viewer
1. Run the following command to initialize the viewer
````bash
blobtools view --remote ~/btk
````
2. Open your browser and go to the URL indicated in the previous command (e.g http://localhost:8001/view/btk)
