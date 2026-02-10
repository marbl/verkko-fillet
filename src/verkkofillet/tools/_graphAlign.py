import sys
import os
import subprocess
import tqdm
import shlex
import pandas as pd
from tqdm import tqdm

def graphIdx(obj,
             graph = "assembly.homopolymer-compressed.gfa",
             threads=1,
             GraphAligner_path="GraphAligner",
             working_directory = "graphAlignment",
             prefix=None, 
             diploid_heuristic_cache = None,
             seeds_mxm_cache = None,
             showOnly = False
             ):
    """\
    Index the graph for graph alignment.

    Parameters
    ----------
    graph
        The path to the graph file. Default is "assembly.homopolymer-compressed.gfa".   
    threads
        The number of threads to use. Default is 1.
    GraphAligner_path
        The path to the GraphAligner executable. Default is "GraphAligner".
    prefix
        The prefix for the output files. Default is None. If None, the graph name is used as the prefix.
    diploid_heuristic_cache
        The path to the diploid heuristic cache file. Default is None. If None, the prefix is used.
    seeds_mxm_cache 
        The path to the seeds mxm cache file. Default is None. If None, the prefix is used.
    showOnly
        If True, the command is printed but not executed. Default is False.
    Returns
    -------
        index files. 
    """

    print(f"Indexing graph: {graph}")
    print(f"But this step will take a while and use a lot of memory in jupyter notebook.")
    print(f"Highly recommend running this step in terminal. If you want to see the command, set showOnly=True")
    print(" ")
    
    graph = os.path.abspath(graph)
    verkkoFilletDir = obj.verkko_fillet_dir
    # print(f"verkkoFilletDir: {verkkoFilletDir}")
    if working_directory == "graphAlignment":
        working_directory = os.path.join(verkkoFilletDir, working_directory)
    else:
        working_directory = os.path.abspath(working_directory)

    # print(f"print pwd : {os.getcwd()}")
    # print(f"working_directory: {working_directory}")
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
        print(f"Folder {working_directory} created.")
    if not os.path.exists(os.path.join(working_directory, 'log')):
        os.makedirs(os.path.join(working_directory, 'log'))
        print(f"Folder {os.path.join(working_directory, 'log')} created.")

    if prefix == None:
        prefix = os.path.basename(graph).replace('.gfa', '')
        diploid_heuristic_cache = f"{working_directory}/{prefix}_diploid-heuristic.index"
        seeds_mxm_cache = f"{working_directory}/{prefix}_seeds-mxm.index"
        
    # check if the graph index exists if so, return
    if os.path.exists(f"{diploid_heuristic_cache}"):
        print(f"{diploid_heuristic_cache} index exists.")
        return
    
    if os.path.exists(f"{seeds_mxm_cache}"):
        print(f"{seeds_mxm_cache} index exists.")
        return
    
    log = os.path.join(verkkoFilletDir, 'log', 'graph_index.log')
    # Proceed only if the index file doesn't exist
    cmd = (
        f"touch {working_directory}/empty.fasta && "
        f"{GraphAligner_path} -t {threads} -g {graph} "
        f"-f {working_directory}/empty.fasta "
        f"-a {working_directory}/empty.gaf "
        f"--diploid-heuristic 21 31 "
        f"--diploid-heuristic-cache {diploid_heuristic_cache} "
        f"--seeds-mxm-cache-prefix {seeds_mxm_cache} "
        f"--bandwidth 15 --seeds-mxm-length 30 --mem-index-no-wavelet-tree "
        f"--seeds-mem-count 10000 > {log} && "
        f"rm {working_directory}/empty.fasta {working_directory}/empty.gaf"
    )

    if showOnly:
        print(f"Command: {cmd}")
        return
    try:
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,  # Captures errors for debugging
            shell=True,
            check=True
        )
        print("Indexing completed successfully!")
        print(f"log file: {log}")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e.stderr.decode().strip()}")


def graphAlign(obj,
               graph = "assembly.homopolymer-compressed.gfa",
               threads=1,
               GraphAligner_path="GraphAligner",
               prefix=None, 
               ReadList = None,
               alignNode = False,
               diploid_heuristic_cache = None,
               seeds_mxm_cache = None,
               working_directory = "graphAlignment",
               showOnly = False):
    """\
    Align ONT reads to the graph.

    Parameters
    ----------
    
    
    Returns
    -------
        alignment files

    """
    # check if the graph index exists
    graph = os.path.abspath(graph)
    if working_directory == "graphAlignment":
        verkkoFilletDir = obj.verkko_fillet_dir
        working_directory = os.path.join(verkkoFilletDir, working_directory)
    else:
        working_directory = os.path.abspath(working_directory)
    print(f"working_directory: {working_directory}")
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
        print(f"Folder {working_directory} created.")
    if not os.path.exists(os.path.join(working_directory, 'log')):
        os.makedirs(os.path.join(working_directory, 'log'))
        print(f"Folder {os.path.join(working_directory, 'log')} created.")

    if prefix == None:
        prefix = os.path.basename(graph).replace('.gfa', '')  
    
    if diploid_heuristic_cache == None:
        diploid_heuristic_cache = f"{working_directory}/{prefix}_diploid-heuristic.index"
        # check if the diploid index exists if not, run the graphIdx function
        if not os.path.exists(diploid_heuristic_cache):
            print("Diploid index does not exist. Running graphIdx function...")
    
    multimapScoreFraction = 0.99
    if alignNode:
        multimapScoreFraction = 0
        print(f"mode : alignNode, multimapScoreFraction = {multimapScoreFraction}")
    
    if seeds_mxm_cache == None:
        seeds_mxm_cache = f"{working_directory}/{prefix}_seeds-mxm.index"
        # check if the seeds-mxm index exists if not, run the graphIdx function
        if not os.path.exists(seeds_mxm_cache):
            print("Seeds-mxm index does not exist. Running graphIdx function...")

    # Construct the graph path within the function
    working_directory = os.path.abspath(working_directory)
    
    # Step 1: Create alignment folder if it doesn't exist
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
        print(f"Folder {working_directory} created.")
    
    # Step 2: Generate ONT read list
    if ReadList==None:
        print(f"ReadList not provided. Generating ReadList from ONT files under 3-align/split/")
        if obj == None:
            print("Please provide the VerkkoFillet object for using the default ReadList generation.")
            return
        ReadList = os.path.join(working_directory, f"{prefix}.ReadList.txt")
        ontReads = os.path.join(obj.verkkoDir, "3-align/split/")
        if not os.path.exists(ontReads):
            print(f"ONT reads directory {ontReads} does not exist.")
            return
        cmd = f"ls {ontReads}*.fasta.gz > {ReadList}"
        print(f"Command: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating ONT read list: {e}")
            return
        
    elif ReadList.endswith((".fa",".fasta",".fasta.gz",".fa.gz",".fastq",".fastq.gz","fq","fq.gz")):
        print(f"single sequence file provided: {ReadList}")
        ReadList = os.path.abspath(ReadList)
        subprocess.run(f"echo {ReadList} > {ReadList}.list", shell=True, check=True)
        ReadList = f"{ReadList}.list"
  
    elif os.path.exists(ReadList):
        ReadList = os.path.abspath(ReadList)

    # Step 3: Align reads
    gaf_path = os.path.join(working_directory, f"{prefix}.gaf")
    
    if os.path.exists(gaf_path):
        print(f"Alignment file {gaf_path} already exists.")
        return

    elif not os.path.exists(gaf_path):
        print(f"Aligning reads to graph: {graph}")
        read_list = pd.read_csv(ReadList, header=None)[0].tolist()
        
        for i in tqdm(range(len(read_list)), desc="Processing reads"):
            read_file = read_list[i]
            gaf_file = os.path.join(working_directory, f"{prefix}_ont{i}.gaf")
            log_file = os.path.join(working_directory, f"{prefix}_ont{i}.log")
        
            # Safely construct the command
            cmd = (
                f"{shlex.quote(GraphAligner_path)} -t {threads} -g {shlex.quote(graph)} "
                f"-f {shlex.quote(read_file)} -a {shlex.quote(gaf_file)} "
                f"--diploid-heuristic 21 31 "
                f"--diploid-heuristic-cache {shlex.quote(diploid_heuristic_cache)} "
                f"--seeds-mxm-cache-prefix {shlex.quote(seeds_mxm_cache)} "
                f"--seeds-mxm-windowsize 5000 "
                f"--seeds-mxm-length 30 --seeds-mem-count 10000 "
                f"--bandwidth 15 --multimap-score-fraction {multimapScoreFraction} "
                f"--precise-clipping 0.85 --min-alignment-score 5000 "
                f"--hpc-collapse-reads --discard-cigar "
                f"--clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 "
                f"--max-trace-count 5 --mem-index-no-wavelet-tree > {shlex.quote(log_file)}"
            )
            if showOnly :
                print(f"Command: {cmd}")
                print(f"This shows the first read only. Please run this command in terminal to align all reads using a loop.")
                break
            else:
                try:
                    subprocess.run(cmd, shell=True, check=True, cwd= working_directory)
                except subprocess.CalledProcessError as e:
                    print(f"Error during alignment of read {i}: {e}")
                    return
        
        # Concatenate GAF files
        concat_cmd = (
            f"cat {os.path.join(working_directory, f'{prefix}_ont*.gaf')} > {gaf_path} && "
            f"rm {prefix}_ont* "
        )
        print(f"Concatenation command: {concat_cmd}")
        
        if not showOnly:
            try:
                subprocess.run(concat_cmd, shell=True, check=True, cwd= working_directory)
            except subprocess.CalledProcessError as e:
                print(f"Error during concatenation: {e}")
                return
            
            print(f"Alignment completed. Final GAF file: {gaf_path}")



def extractNodeSeq(node,
                   graph = 'assembly.homopolymer-compressed.gfa',
                   working_directory = "graphAlignment",
                   outPrefix = None):
    """\
    Extract the sequence of a node from the graph.

    Parameters
    ----------
    node
        The node to be extracted.
    graph
        The path to the graph file. Default is "assembly.homopolymer-compressed.gfa".
    working_directory
        The directory to store the extracted sequence. Default is "graphAlignment".
    outPrefix
        The prefix for the output file. Default is None. If None, the node name is used as the prefix.

    Returns
    -------
        sequence file
    """
    if outPrefix is None:
        outPrefix = node
        print(f"Output prefix not provided. Using {outPrefix} as output prefix")

    print(f"Extracting sequence for {node}...")
    working_directory = os.path.abspath(working_directory)
    if not os.path.exists(working_directory):
        print(f"Creating directory {working_directory}")
        os.makedirs(working_directory)

    graph = os.path.abspath(graph)

    cmd = f"awk -v node={node} '$1 == \"S\" && $2 == node {{print \">\"node\"\\n\"$3}}' {graph} > {outPrefix}.fa"

    try:
        subprocess.run(cmd, shell=True, check=True, cwd= working_directory)
        print(f"Sequence extracted successfully. Output file: {working_directory}/{outPrefix}.fa")
    except subprocess.CalledProcessError as e:
        print(f"Error during extracting sequence of {node}: {e}")