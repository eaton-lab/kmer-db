# km-unity
A community sourced database of genome characteristics (genome size and heterozygosity estimates) extracted from published Illumina data. The goal of this project is to:


1. Provide a simple tool for downloading Illumina data from NCBI and extracting kmer-based statistics in a unified workflow.
2. The tool is easy to use, cleans up after itself, and is easily distributed among users and on HPC.
3. Data are publicly available, continually updated and community maintained through GitHub pull requests.
4. Log files store details of all analyses.
5. Samples can be selected by user or auto-populated by querying NCBI for unpopulated results.


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



| Organism  | Taxid  |  Biosample  |   Run   |  Bases_Gb  |   Estimate_Genome_Size   |  Estimate_Heterozygosity  |
|   ---     |   ---  |     ---     |   ---   |    ---     |          ---             |          ---              | 
|   Trachypithecus laotum         |   465718 |   SRS5312446    |  SRR10028098   |  116   |   ...    |   ...   |
|   Trachypithecus poliocephalus  |   34886  |   SRS3420064    |  SRR7345456    |  121   |   ...    |   ...   |



| Biosample  | Organism |  Tax_id     |  Size_Gb   |  SRR   |   genome_size   |  heterozygosity  |
| ---        |  ---     |   ---       |   ---      |  ---   |   ---   |  ---   |
| SRS5312446 | Trachypithecus laotum 	| 465718 	|  116 	|  SRR10028098 	|  0.513780  |
| SRS3420064 | Eumetopias jubatus 	  | 34886 	|  121 	|  SRR7345456 	|  0.295698  |
| SRS3734713 | Trachypithecus poliocephalus     |	465719 | 150 | SRR7778912 |  0.509932  |
| SRS2059369 | Ammotragus lervia      | 9899 	  |   76   |	SRR5438049  |	0.221502 |
...
