# BTKCourse
## BlobToolKit (v3.3.0) instaltion 

````bash
pip install blobtoolkit
````

The blobtools viewc ommand requires firefox or a chromium-based browser to start the interactive viewer or to generate plots from the command line, these can be installed with:

````bash
conda install -c conda-forge firefox geckodriver
````

## Demo on BTK
```mermaid
graph TB
subgraph CASE 2
    er(Raw Reads) -->|extracte target reads| aa(Newdenovo assembly)
    aa -- <b>Blobtoolkit --> t(Clean Assembly)
   
  end
  subgraph CASE 1
    tt(Raw Reads) --> a
    a[denovo assembly] --> rr(<i>Blobtoolkit)
    rr --> b(Contaminated assembly!)
   

    rr --> zz(Clean Assembly!)
    b -->|filter| EE(target group contigs)
  end
  
        classDef green fill:#93FF33,stroke:#333,stroke-width:2px
        classDef blue fill:#00FA9A,stroke:#333,stroke-width:4px
       
        class g,a,h green
        class b,c,d,e,f blue
 ```      


## Create a project directory

````bash
mkdir btk
````
All following commands will executed inside this folder

## Fetch Databases

### 1. Fetch the NCBI Taxdump
````bash
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
````
### 2. Fetch the nt database

````bash
mkdir -p nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P nt/ && \
        for file in nt/*.tar.gz; \
            do tar xf $file -C nt && rm $file; \
        done
````
### 3. Fetch and format the UniProt database

````bash
mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')
cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo -e "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
cd -
````


## Create blobtoolkit project/database



## Input file requirements

metadata file in yaml format, to describe the project e. g. B_cinera.yaml, Can be downloaded from HERE.
one assembly file, e.g. assembly.fasta, Can be downloaded from HERE.
one (or more) coverage file(s) e.g. mapping_1.bam, Can be downloaded from HERE.
one (or more) hits file(s), e.g. blastn.out and diamond.blastx.out, Can be downloaded from HERE.


### 1. Create a metadata file

````bash
vim B_cinera.yaml
````
Write the following information inside the yaml file. 

````bash
assembly:
  alias: B_cinera_112
  record_type: contig
taxon:
  name: Botrytis cinerea
````

### 2. Create a hit file(s)

#### 2.1. The blastn hit file 

````bash
blastn -db ./nt/nt \
-query assembly.fasta \
-outfmt "6 qseqid staxids bitscore std" \
-max_target_seqs 10 -max_hsps 1 \
-evalue 1e-25 -num_threads 30 \
-out blastn.out
````

#### 2.2. The deamond hit file 

````bash
diamond blastx --query assembly.fasta \
--db ./uniprot/reference_proteomes.dmnd \
--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
--sensitive --max-target-seqs 1 \
--evalue 1e-25 \
--threads 30 > diamond.blastx.out
````

### 3. Create a coverage file(s)

mapping file using minimap


````bash
minimap2 -ax sr -t 30 assembly.fasta \
trimmed.R1.fastq trimmed.R2.fastq \
| samtools sort -@30 -O BAM -o coverage.bam -
````

### 3. Create a BUSCO summary file

Run busco on the genome assemlby 

````bash
busco -i assembly.fasta \
-l helotiales_odb10 -o botrytis -m genome --cpu 30
````


### 1. Create a metadata file

````bash
blobtools create --fasta Assembly.fasta \
--meta Assembly.yaml  --taxid 75913 \
--taxdump /home1/software/blobtoolkit/taxdump ~/btk
````

### 2. Adding hits
3. Add hits with blobtools command:
````bash
blobtools add --hits asm.ncbi.blastn.out \
--hits asm.diamond.blastx.out \
--taxrule bestsumorder \
--taxdump /home1/software/blobtoolkit/taxdump Deconta_asm/
````
### 3. Adding Coverage

2. add coverage using blobtools command:
````bash
blobtools add --cov asm.bam --threads 30 ~/btk
```` 
### 4. Adding buscos
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
