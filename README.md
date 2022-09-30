# BTKCourse
## BlobToolKit (v3.3.0) instaltion 

````bash
pip install blobtoolkit
````

The blobtools viewcommand requires firefox or a chromium-based browser to start the interactive viewer or to generate plots from the command line, these can be installed with:

````bash
conda install -c conda-forge firefox geckodriver
````

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
## Databases
# 1. Fetch the NCBI Taxdump
````bash
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
````
# 2. Fetch the nt database

````bash
mkdir -p nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P nt/ && \
        for file in nt/*.tar.gz; \
            do tar xf $file -C nt && rm $file; \
        done
````
# 3. Fetch and format the UniProt database

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

## Adding data to a dataset (Create blobtoolkit project)
### 1. Create a metadata file
1. Create a project directory
````bash
mkdir btk
````
2. Create yaml file to describe the project:

````bash
vim Assembly.yaml
````
Write the following information inside the yaml file. You download the Assembly.yaml from HERE.

````bash
assembly:
  alias: B_cinera_112
  record_type: contig
taxon:
  name: Botrytis cinerea
````

3. Run the following command to creat btk directory project

````bash
blobtools create --fasta Assembly.fasta \
--meta Assembly.yaml  --taxid 75913 \
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
minimap2 -ax sr -t 30 /sanhome2/Comparative_genomics/Botrytis/Assembly/denovoCLC/Contig_112-name.trim.asco.fa \
../kraken/112-name.trim.asco.R1.fastq ../kraken/112-name.trim.asco.R2.fastq \
| samtools sort -@30 -O BAM -o asm.DRR008460.bam -
````
2. add coverage using blobtools command:
````bash
blobtools add --cov asm.bam --threads 30 ~/btk
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
