# kmunity (Under development!)
A community sourced database of genome characteristics (genome size and heterozygosity estimates) extracted from published Illumina data. The goal of this project is to provide a simple tool to:


#### New TODO: Move to GCP Bigquery
- https://www.youtube.com/watch?v=PoseTXa8LU8&list=PLH-TjWpFfWrt5MNqU7Jvsk73QefO3ADwD&index=3  
- Data transfer with buckets  
- Avoid slow and expensive need for data egress (downloads)  
- ...


1. Query NCBI Run accession IDs (SRR) to fetch data and metadata. 
2. Downloading fastq data from NCBI (in a transparent way that does not leave behind zombies like sra-tools.)
3. Calculate kmer statistics from data with kmerfreq, gce, etc. 
4. Organize results into a public repository on GitHub.
5. Record detailed log files for all contributed analyses.


### Repository structure
```bash
├── kmunity
│   └── kmunity
│       ├── __init__.py
│       └── Kmunity.py
├── mammals
│   ├── database.csv
|   └── logfiles
|       ├── SRRXXXYYY.log
|       └── ...
├── birds
│   ├── database.csv
|   └── logfiles
|       ├── SRRXXXYYY.log
|       └── ...
├── other
│   ├── database.csv
|   └── logfiles
|       ├── SRRXXXYYY.log
|       └── ...
└── plants
    ├── database.csv
    └── logfiles
        ├── SRRXXXYYY.log
        └── ...
```


### database files

| Organism  | Taxid  |  Biosample  |   Run   |  Bases_Gb  |  Coverage |  Genome_Size   |  Heterozygosity  |
|   ---     |   ---  |     ---     |   ---   |    ---     |          ---             |          ---     |   ---         | 
|   Trachypithecus laotum         |   465718 |   SRS5312446    |  SRR10028098   |  116   |   ...    |   ...   |   ... |
|   Trachypithecus poliocephalus  |   34886  |   SRS3420064    |  SRR7345456    |  121   |   ...    |   ...   |  ... |


### Installation (linux only)
```bash
# install kmunity
git clone https://github.com/{username}/kmunity
cd kmunity/
pip install .

# configure sra-tools >=v.2.10.5
kmunity --config
```


### Contributing to database
```bash
# run tool selecting desired database to contribute to, scratch space, and repo location
kmunity -d mammals  -w /scratch/tmp  -r ./kmunity

# or select a specific run ID and kmunity will auto-select the appropriate database
kmunity -s SRR10028098 -w /scratch/tmp  -r ./kmunity

# cd into dir, diff to see new results, push changes, and make pull request on GitHub.
cd kmunity/
git diff 
git add *
git commit -m "added stats for SRR10028098"
git push
```
