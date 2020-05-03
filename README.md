# kmunity (Under development!)
A community sourced database of genome characteristics (genome size and heterozygosity estimates) extracted from published Illumina data. The goal of this project is to:

1. Provide a simple tool for downloading Illumina data from NCBI and extracting kmer-based statistics in a unified workflow.
2. The tool is easy to use, cleans up after itself, and is easily distributed among users and on HPC.
3. Data are publicly available, continually updated and community maintained through GitHub pull requests.
4. Log files store details of all analyses.
5. Samples can be selected by user or auto-populated by querying NCBI for unpopulated results.


### Repository structure
```bash
├── src
│   ├── kmerhet
│       ├── __init__.py
│       └── KmerHet.py
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
├── plants
│   ├── database.csv
|   └── logfiles
|       ├── SRRXXXYYY.log
|       └── ...
```


### database files

| Organism  | Taxid  |  Biosample  |   Run   |  Bases_Gb  |   Estimate_Genome_Size   |  Estimate_Heterozygosity  |
|   ---     |   ---  |     ---     |   ---   |    ---     |          ---             |          ---              | 
|   Trachypithecus laotum         |   465718 |   SRS5312446    |  SRR10028098   |  116   |   ...    |   ...   |
|   Trachypithecus poliocephalus  |   34886  |   SRS3420064    |  SRR7345456    |  121   |   ...    |   ...   |


### Contributing to database
```bash
# fork and clone the repo
git clone https://github.com/{username}/kmunity

# install kmunity
conda install kmunity -c conda-forge

# run tool selecting desired database to contribute to, scratch space, and repo location for results.
# multiple jobs can be distributed on HPC at same time if you have sufficient scratch space.
kmunity -db mammals  -o /scratch/tmp  --repo ./kmunity

# cd into dir, check new results, push, and make pull request to origin on GitHub
cd kmunity
git diff 
git push
```
