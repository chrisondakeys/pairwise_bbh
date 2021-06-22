# pairwise_bbh
Finds orthologous genes given a query proteome and a pool of proteomes. 
Say file a.faa is a query proteome and dir a directory containing a number of proteomes, pairwise_bbh will run a bidirectional best hit using a pairwise blastp (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) search between every sequence in a.faa and each proteome present in dir. Once done, a custom output directory will be populated with a set of ordered multiFASTAs for each of the input proteomes (query included). A semicolon-delimited summary table is printed to stdout, so you may want to run this script as:  

```python3 pairwise_bbh.py -q a.faa -d dir -o Output > table.csv```  

Note that the input query (a.faa in this example) will remain untouched and a shrinked copy moved to the output directory.  

## Dependencies 
- This version of pairwise_bbh runs on Linux and MAC OSX, and it uses native shell commands throughout the computation. 
- The only external dependency is the blast+ suite, or at least ```makeblastdb``` and ```blastp```. 
- Read and write privileges.  

### Usage  
```Usage: phyloprep.py [OPTIONS] -q <file> -d <path> -o <str>  
Mandatory args:  
	-q | --query <file>	protein (multi)FASTA to be sampled. Will blast each sequence against the database  
	-d | --database <path>	directory containing proteomes against which orthologs will be searched  
	-o | --output <str>	output directory where results will be stored in. Must not be an existing path  
Optional args:  
	-e | --evalue <float>	cutoff value for dropping hits. default = 1e-30  
	-t | --threads <int>	number of threads for blastp. default = 1  
  ```
  
  ## Example  
  
  After downloading and uncompressing the material at the link that follows, you may change directory into the master folder and run pairwise_bbh.py as follows:  
  
  
```python3 pairwise_bbh.py -q Pseudoalteromonas_arctica_shrinked.faa -d Database -o Output > table.csv``` 

This will use all the sequences present in Pseudoalteromonas_arctica_shrinked.faa as BBH queries against all the genomes present in the Database directory. As a result, you should find the Output directory populated with as many proteomes there were in the Database directory, including the query. Each file should contain sequences that are orthologs to the initial query. By writing stdout to a table.csv, you will also store a semicolon-delimited file containing the protein identifiers that were found for each organism.  

  
[pairwise_bbh-master.zip](https://github.com/chrisondakeys/pairwise_bbh/files/6695953/pairwise_bbh-master.zip)

