"""
===========================================================================================================================
* AUTHOR:		Christopher Riccardi, PNRA-funded postgrad reasercher at Computational Biology Dept., University of Firenze
* PROJECT:		Identification of transcriptional regulatory networks in Antarctic bacteria as a proxy for global warming 
				effects on microbial life in the Polar Oceans
* PURPOSE:		Find genes that are orthologous to a query, given a pool of other species
===========================================================================================================================
"""
import multiprocessing
import random
import string
import subprocess
import sys
import os

global CHARS

CHARS = 16
TOTAL_BLAST_SEARCHES = 0

class CommandLine(object):
    """
    Arguments:
        -q | --query
        -d | --database
        -o | --output
        -e | --evalue
        -t | --threads
    """
    def __init__(self, argv):
        self.argv = argv
        self.argc = len(argv)
        self.eval_cutoff = 1e-30
        self.input_query = ""
        self.input_database = ""
        self.output_dir = ""
        self.threads = 1
        self.max_threads = multiprocessing.cpu_count()
        global EVALUE_THRESH
        global THREADS
        EVALUE_THRESH = self.eval_cutoff
        THREADS = self.threads
    def print_usage(self):
        q_arg = "\t-q | --query <file>\tprotein (multi)FASTA to be sampled. Will blast each sequence against the database\n"
        i_arg = "\t-d | --database <path>\tdirectory containing proteomes against which orthologs will be searched\n"
        o_arg = "\t-o | --output <str>\toutput directory where results will be stored in. Must not be an existing path\n"
        e_arg = "\t-e | --evalue <float>\tcutoff value for dropping hits. default = 1e-10\n"
        t_arg = "\t-t | --threads <int>\tnumber of threads for blastp. default = 1\n"
        
        sys.stderr.write("Usage: " + self.argv[0] + " [OPTIONS] -q <file> -d <path> -o <str>\n")
        sys.stderr.write("Mandatory args:\n")
        sys.stderr.write(q_arg)
        sys.stderr.write(i_arg)
        sys.stderr.write(o_arg)
        
        sys.stderr.write("Optional args:\n")
        sys.stderr.write(e_arg)
        sys.stderr.write(t_arg)
        sys.stderr.write("\n")
    def parse(self):
        if self.argc < 7 or self.argc > 12:
            self.print_usage()
            sys.exit(1)
        for i in range(1, self.argc - 1):
            if self.argv[i] == "-q" or self.argv[i] == "--query":
                if(os.path.isfile(self.argv[i + 1])):
                    self.input_query = os.path.abspath(self.argv[i + 1])
                    i += 2
                else:
                    self.print_usage()
                    sys.exit(1)
                       
            elif self.argv[i] == "-d" or self.argv[i] == "--database":
                if(os.path.isdir(self.argv[i + 1])):
                    self.input_database = self.argv[i + 1]
                    i += 2
                else:
                    self.print_usage()
            elif self.argv[i] == "-o" or self.argv[i] == "--output":
                if(os.path.isdir(self.argv[i + 1])): 
                    ## throw path exists error
                    print("[ERROR] path already exist. Please provide another name", file = sys.stderr)
                    self.print_usage()
                    sys.exit(1)
                    
                else:
                    self.output_dir = os.path.abspath(self.argv[i + 1])
                    global OUTPUT_DIR
                    OUTPUT_DIR = self.output_dir
                    i += 2
            elif self.argv[i] == "-e" or self.argv[i] == "--evalue":
                try:
                    self.eval_cutoff = float(self.argv[i + 1])
                    global EVALUE_THRESH
                    EVALUE_THRESH = self.eval_cutoff
                    i += 2
                except:
                    self.print_usage()
                    sys.exit(1)
                    ## throw float not typed error
            elif self.argv[i] == "-t" or self.argv[i] == "--threads":
                try:
                    self.threads = int(self.argv[i + 1])
                    global THREADS
                    THREADS = self.threads
                    i += 2
                except:
                    ## throw integer not typed error
                    print("[ERROR] threads must be a valid integer", file = sys.stderr)
                    self.print_usage()
                    sys.exit(1)
                if self.threads < 1 or self.threads > self.max_threads:
                    print(F"[ERROR] threads must be an integer between 1 and {self.max_threads}", file = sys.stderr)
                    self.print_usage()
                    sys.exit(1)
                    
class Matrix(object):
    def __init__(self):
        #rows, cols = i, j
        self.cols = []
        self.rows = []
        self.product = []
    def cross(self):
        for j in self.rows:
            for i in self.cols:
                if j != i:
                    self.product.append([j, i])
    
class Index(object):
    def __init__(self, path):
        self.path = path
        self.elems = []
        self.sizes = []
        self.pairs = []
    def fill(self):
        for listdir in os.listdir(self.path):
            if(listdir.endswith(".faa")):
                #proteome = os.path.join(self.path, listdir)
                proteome = listdir
                self.elems.append(os.path.abspath(proteome))
                prot_size = 0
                with open(proteome, "r") as f:
                    for line in f:
                        if(line.startswith('>')):
                            prot_size += 1
                self.sizes.append(prot_size)
    def pop_random(self, lst):
        idx = random.randrange(0, len(lst))
        return lst.pop(idx)
    def random_couples(self):
        elems = self.elems
        if len(elems) % 2 != 0 and len(elems) > 2:
            elems.pop(len(elems) - 1)
        while(elems):
            rand1 = self.pop_random(elems)
            rand2 = self.pop_random(elems)
            pair = [rand1, rand2]
            self.pairs.append(pair)
    def square_matrix(self):
        M = Matrix()
        M.rows = self.elems
        M.cols = self.elems
        M.cross()
        return M.product
    def adjacency_list(self, query):
        M = Matrix()
        M.rows.append(query)
        M.cols = self.elems
        M.cross()       
        return M.product

class FASTA(object):
    def __init__(self, path):
        self.path = path
        self.size = 0
        
        self.sequences = []
        buf = ""
        printing = False
        with open(self.path, "r") as m:
            for line in m:
                if line.startswith(">"):
                    self.size += 1
                    if printing:
                        self.sequences.append(buf)
                        printing = False
                        buf = ""
                        buf += line
                    else:
                        buf += line
                else:
                    buf += line
                    printing = True
            self.sequences.append(buf)
    def seq_by_ID(self, string):
        for j in self.sequences:
            tmp = j.split(">")[1].split()[0].rstrip()
            if string == tmp:
                return j
        return ""
                    
            
## functions
def pairwise_complete_BBH(prot1, prot2, shrinked):
    ## takes two FASTA objects as input
    ## make blast database
    if not os.stat(prot1.path).st_size or not os.stat(prot2.path).st_size:
        sys.stderr.write("[DONE] Sorry, your sequences do not contain any orthologs\n", file = sys.stderr)
        return 1

    ## temporary file
    TMP_FILE = "tmp"+"".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(CHARS))  

    copy_for_shrink = []    ## store indices
    for j in prot1.sequences:
        f = open(TMP_FILE, "w")
        f.write(j)
        f.close()
        blast_cmd1 = "blastp -query " + TMP_FILE + " -db " + prot2.path + " -num_threads " + str(THREADS) + " -evalue " + str(EVALUE_THRESH) + "| grep -A2 'Sequences producing significant alignments:' | tail -n 1" 
        stdout_seq2 = subprocess.check_output(blast_cmd1, shell=True)
        stdout_seq2 = stdout_seq2.decode("utf-8")

        if len(stdout_seq2) < 1:
            ## blastp did not produce any result
            os.remove(TMP_FILE)
            copy_for_shrink.append(prot1.sequences.index(j))
        else:
            ## blastp found a hit
            initial_identifier = j.split(">")[1].split()[0].rstrip()
            protein_id_seq2 = stdout_seq2.split()[0].rstrip()
            os.remove(TMP_FILE)
            ## once found, retrieve the entire sequence and write it to file
            ## in order for it to be blastp'd against seq1, this is the very
            ## core of BBH
  
            hit_sequence = prot2.seq_by_ID(protein_id_seq2)
            if hit_sequence == "":
                sys.stderr.write("Sorry, something went wrong with " + protein_id_seq2 + "\n")
                return 1
            f = open(TMP_FILE, "w")
            f.write(hit_sequence)
            f.close()
            ## now the actual blastp
            blast_cmd2 = "blastp -query " + TMP_FILE + " -db " + prot1.path + " -num_threads " + str(THREADS) + " -evalue " + str(EVALUE_THRESH) + "| grep -A2 'Sequences producing significant alignments:' | tail -n 1" 
            stdout_seq1 = subprocess.check_output(blast_cmd2, shell=True)
            stdout_seq1 = stdout_seq1.decode("utf-8")
            
            if len(stdout_seq1) < 1:
                ## blastp did not produce any result
                os.remove(TMP_FILE)
                copy_for_shrink.append(prot1.sequences.index(j))
            else:
                ## blastp found a hit
                protein_id_seq1 = stdout_seq1.split()[0].rstrip()
                os.remove(TMP_FILE)
                ## compare to initial sequence identifier
                if protein_id_seq1 != initial_identifier:
                    copy_for_shrink.append(prot1.sequences.index(j))
    ## shrink data
    if len(copy_for_shrink) == 0:
        return 0
    #os.remove(prot1.path)
    print(F"Shrinked {len(copy_for_shrink)} sequences", file=sys.stderr)
    f = open(shrinked, "w")
    for i in range(len(prot1.sequences)):
        if i not in copy_for_shrink:
            f.write(prot1.sequences[i])
    f.close()
    return 0
    

def finalize_results(query, index_obj):
    rmdb_cmd = "rm *.pdb *.phr *.pin *.pot *.psq *.ptf *.pto 2> /dev/null"
    os.system(rmdb_cmd)
    query_sequence = FASTA(query).sequences
    mkdb_cmd = "makeblastdb -dbtype prot -in " + query
    subprocess.check_output(mkdb_cmd, shell=True)
    
    for elem in index_obj.elems:
        mkdb_cmd = "makeblastdb -dbtype prot -in " + elem
        subprocess.check_output(mkdb_cmd, shell=True)
    
    print("Finalizing results", file = sys.stderr)
    ## temporary file
    TMP_FILE = "tmp"+"".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(CHARS))  

    headers = [query]
    df = []
    for j in query_sequence:
        df_V1 = [j.split(">")[1].split()[0].rstrip()]
        f = open(TMP_FILE, "w")
        f.write(j)
        f.close()
        
        for i in index_obj.elems:
            if j != i:
                blast_cmd1 = "blastp -query " + TMP_FILE + " -db " + i + " -num_threads " + str(THREADS) + " -evalue " + str(EVALUE_THRESH) + "| grep -A2 'Sequences producing significant alignments:' | tail -n 1" 
                stdout_seq2 = subprocess.check_output(blast_cmd1, shell=True)
                stdout_seq2 = stdout_seq2.decode("utf-8")
        
                if len(stdout_seq2) < 1:
                    ## blastp did not produce any result
                    df_V1 = []
                    break
                else:
                    df_V1.append(stdout_seq2.split()[0].rstrip())
                    if i not in headers:
                        headers.append(i)
        os.remove(TMP_FILE)
        if df_V1:
            df.append(df_V1)
            df_V1 = []

    os.system(rmdb_cmd)
    table_headers = []
    i = 0
    for elem in headers:
        f = open(os.path.splitext(os.path.basename(elem))[0]+".orth", "w")
        table_headers.append(os.path.splitext(os.path.basename(elem))[0])
        j = FASTA(elem)
        for sequence in df:
            f.write(j.seq_by_ID(sequence[i]))
        f.close()
        os.remove(elem)
        i += 1
        
    print(';'.join(table_headers))
    for vector in df:
        print(';'.join(vector))        
        
    path = query.replace(os.path.abspath(query).split("/")[-1], "")
    rmdb_cmd = "rm "
    rmdb_cmd += path + "/*.pdb " + path + "/*.phr " + path + "/*.pin " + path + "/*.pot " + path + "/*.psq " + path + "/*.ptf " + path + "/*.pto 2> /dev/null"
    os.system(rmdb_cmd)

## ******* PROCEDURAL PROGRAMMING ******* ##

## parse input                
args = CommandLine(sys.argv)
args.parse()

## log basic info
print(F"{sys.argv[0]} was called as:\n", file = sys.stderr)
print(" ".join(args.argv), file = sys.stderr)
print(F"\nUsing {THREADS} threads", file = sys.stderr)
print(F"E value threshold set to {EVALUE_THRESH}", file = sys.stderr)

## preliminary operations
os.system("mkdir " + args.output_dir)
ls_cmd = "ls " + args.input_database + " | wc -l 2> /dev/null"
if int(subprocess.check_output(ls_cmd, shell=True)==0):
    sys.exit(1)
cp_cmd = "cp " + args.input_database + "/* " + args.output_dir + " 2> /dev/null"
ls_cmd = "ls " + args.output_dir + " | wc -l 2> /dev/null"
subprocess.check_output(cp_cmd, shell=True)
n_files = subprocess.check_output(ls_cmd, shell=True)
print(f"Copied {int(n_files)} files", file = sys.stderr)

## create blast db
mkdb_cmd = "makeblastdb -dbtype prot -in " + args.input_query
subprocess.check_output(mkdb_cmd, shell=True)

## change working directory
os.chdir(args.output_dir)

for listdir in os.listdir(args.output_dir):
    mkdb_cmd = "makeblastdb -dbtype prot -in " + listdir
    subprocess.check_output(mkdb_cmd, shell=True)


## index elems and evaluate / shrink
index = Index(os.getcwd())
index.fill()
index.pairs = index.adjacency_list(args.input_query)
shrinked = os.path.splitext(os.path.basename(args.input_query))[0]+".shrinked"
shrink_cmd = "cp " + args.input_query + " " + shrinked
os.system(shrink_cmd)
for couples in index.pairs:
    print(f"Pairwise shrinkage on {couples}", file = sys.stderr)
    seq1 = FASTA(couples[0])
    seq2 = FASTA(couples[1])
    pairwise_complete_BBH(seq1, seq2, shrinked)

    
finalize_results(shrinked, index)
path = args.input_query.replace(os.path.abspath(args.input_query).split("/")[-1], "")
rmdb_cmd = "rm "
rmdb_cmd += path + "/*.pdb " + path + "/*.phr " + path + "/*.pin " + path + "/*.pot " + path + "/*.psq " + path + "/*.ptf " + path + "/*.pto 2> /dev/null"
os.system(rmdb_cmd)