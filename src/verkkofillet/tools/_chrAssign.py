import sys
import shlex
import subprocess
import os
import re
import shutil
from IPython.display import Image, display
from datetime import datetime
from .._run_shell import run_shell

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def chrAssign(obj, ref, working_directory="chromosome_assignment", fasta="assembly.fasta", chr_name="chr", idx=99, showOnly=False, min_Length = 5000000, force = False):
    """\
    Run the script to align the assembly to the given reference using mashmap and obtain the chromosome assignment results.

    Parameters:
    -----------
    obj (verko-fillet object):
        An object that contains a .stats attribute, which should be a pandas DataFrame.
    ref (str) :
        Existing reference
    fasta (str):
        verkko assembly. [default: `assembly.fasta`]
    working_directory (str):
        output directory [default : `./stats/`]
    chr_name (str):
        prefix of the chromosome name in the previous reference. [default : "chr"]
    idx (int):
        Identity threshold to filter mashmap result [defualt : 99]
    showOnly (bool): 
        If set to True, the script will not be executed; it will only display the intended operations. [default : FALSE]
    min_Length (int):
        Minimum length of the contigs to be considered for chromosome assignment. [default : 5000000]
    force (bool):
        If set to True, the script will overwrite the existing files. [default : FALSE]

    Return:
    -----------
    {working_directory}/assembly.mashmap.out
    {working_directory}/assembly.mashmap.out.filtered.out
    {working_directory}/chr_completeness_max_hap1
    {working_directory}/chr_completeness_max_hap2
    {working_directory}/translation_hap1
    {working_directory}/translation_hap2
    """
    # Ensure absolute paths
    working_dir = os.path.abspath(working_directory)
    script = os.path.abspath(os.path.join(script_path, "getChrNames.sh"))  # Assuming script_path is defined elsewhere
    fasta = os.path.abspath(fasta)
    ref = os.path.abspath(ref)
    cwd = os.getcwd()
    hcp_ref = f"{os.path.splitext(ref)[0]}.hpc.fasta"

    if not os.path.exists(ref):
        print(f"Reference file not found: {ref}")
        return
        
    # remove all results if force is True
    if force:
        print("Force is set to True. Removing all existing results.")
        print(f"Rerun the job to get the new results.")
        if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
    else:
        if os.path.exists(f"{working_directory}/assembly.mashmap.out"):
            print(f"The mashmap file {working_directory}/assembly.mashmap.out already exists.")
            print(f"If you want to re-run this job, please delete {working_directory}/assembly.mashmap.out or set force=True")
            print(f"Reusing the existing results for this time.")
            shutil.copy(f"{working_directory}/assembly.mashmap.out", f"{cwd}/assembly.mashmap.out")
         
    output_files = [
        "translation_hap1",
        "translation_hap2",
        "chr_completeness_max_hap1",
        "chr_completeness_max_hap2",
        "assembly.mashmap.out",
        "assembly.homopolymer-compressed.chr.csv",
    ]

    # Check if all output files already exist
    if all(os.path.exists(os.path.join(working_dir, file)) for file in output_files):
        print("All output files already exist. Skipping chromosome assignment.")
        return

    cleanupFiles = [
        "/8-hicPipeline/unitigs.hpc.mashmap.out",
        "unitigs.hpc.mashmap.out",
        hcp_ref,
    ] + output_files

    for file in cleanupFiles:
        if os.path.exists(file):
            os.remove(file)

    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    # Check if the working directory exists
    if not os.path.exists(working_dir):
        print(f"Creating the working directory: {working_dir}")
        os.mkdir(working_dir)

    # Construct the shell command
    cmd = f"bash {shlex.quote(script)} {shlex.quote(ref)} {shlex.quote(str(idx))} {shlex.quote(fasta)} {shlex.quote(chr_name)} {shlex.quote(str(min_Length))}"
    
    run_shell(cmd, wkDir=os.getcwd(), functionName="chrAssign", longLog=False, showOnly=showOnly)

    for output in output_files:
        shutil.move(output, f"{working_dir}/{output}")

    # cleanup
    for file in cleanupFiles:
        if os.path.exists(file):
            os.remove(file)



def convertRefName(fasta, map_file, out_fasta=None, showOnly=False):
    """
    Replace the name in the given FASTA file.

    Parameters
    ----------
    fasta
        FASTA file in which the contig name is to be replaced
    map_file
        A two-column file, delimited by tabs, containing the old and new contig names.
    showOnly
        If set to True, the script will not be executed; it will only display the intended operations. [default : False]

    Returns
    -------
        Output fasta file with the new contig names.
    """
    # Default out_fasta if not provided
    ref_fasta = os.path.abspath(fasta)
    map_file = os.path.abspath(map_file)
    working_dir = os.path.dirname(ref_fasta)
    
    if out_fasta is None:
        # Extract the base name of the file and directory
        basename_fasta = os.path.basename(ref_fasta)
        dir_fasta = os.path.dirname(ref_fasta)
        
        # Always remove the last extension (e.g., .gz, .fa, .fasta)
        if basename_fasta.endswith(".gz"):
            basename_fasta = os.path.splitext(basename_fasta)[0]  # Remove .gz
        
        basename_fasta = os.path.splitext(basename_fasta)[0]  # Remove the actual file extension
        out_fasta = os.path.join(dir_fasta, f"{basename_fasta}.rename.fa.gz")

    # Check if the output file already exists
    if os.path.exists(out_fasta):
        print(f"The renamed reference fasta already exists: {out_fasta}")
        return
    
    # Construct the awk command to replace headers
    script = os.path.abspath(os.path.join(script_path, "changeChrName.sh"))
    cmd=f"bash {shlex.quote(script)} {shlex.quote(map_file)} {shlex.quote(ref_fasta)} {shlex.quote(out_fasta)}"

    if showOnly:
        # If showOnly is True, just display the command instead of executing it
        print(f"Command to be executed:\n{cmd}")
    else:
        # Run the shell command to perform the operation
        run_shell(cmd, wkDir=working_dir, functionName="convertRefName", longLog=False, showOnly=showOnly)


def showPairwiseAlign(obj, 
                      size: str ="large", 
                      working_directory: str ="chromosome_assignment",
                      mashmap_out: str ="chromosome_assignment/assembly.mashmap.out", 
                      prefix: str ="refAlign", 
                      idx: float =0.99, 
                      minLen: int =50000, 
                      showOnly: bool =False):
    """
    Generate a dot plot from the mashmap output.

    Parameters
    ----------
    obj
        An object that contains a .stats attribute, which should be a pandas DataFrame.
    size
        Size of the image.
    mashmap_out
        Path to the mashmap output file.
    prefix
        Prefix for the output files.
    idx
        Identity threshold for filtering alignments.
    minLen
        Minimum length of alignments to be considered.
    showOnly
        If set to True, the script will not be executed; it will only display the intended operations.
    """
    # Ensure paths are absolute
    working_dir = os.path.abspath(working_directory)
    script = os.path.abspath(os.path.join(script_path, "generateDotPlot"))
    log_file = os.path.join(working_dir, "logs", "showPairwiseAlign.log")
    mashmap_out = os.path.abspath(mashmap_out)
    
    # Check if gnuplot is available
    gnuplot_path = shutil.which("gnuplot")
    if not gnuplot_path:
        print(f"Command 'gnuplot' is not available.")
        return

    # Filtering command
    cmd1 = (
        f"awk -F'\t' '{{ split($13, arr, \":\"); "
        f"if ((arr[3] > {idx}) && ($4 - $3 > {minLen})) print }}' "
        f"{shlex.quote(mashmap_out)} > {shlex.quote(mashmap_out)}.filtered.out"
    )
    run_shell(cmd1, wkDir=working_dir, functionName="showPairwiseAlign_1", longLog=False, showOnly=showOnly)

    # Generate plot command
    cmd2 = f"perl {shlex.quote(script)} png {shlex.quote(size)} {shlex.quote(mashmap_out)}.filtered.out"
    run_shell(cmd2, wkDir=working_dir, functionName="showPairwiseAlign_2", longLog=False, showOnly=showOnly)

    # Rename output files
    output_files = ['out.fplot', 'out.rplot', 'out.gp', 'out.png']
    for file in output_files:
        old_path = os.path.join(working_dir, file)
        new_path = os.path.join(working_dir, f"{prefix}.{file.split('.')[1]}")
        if os.path.exists(old_path):
            os.rename(old_path, new_path)
    
    # Display the PNG image
    image_path = os.path.join(working_dir, f"{prefix}.png")
    if os.path.exists(image_path):
        display(Image(filename=image_path, width=500))
    else:
        print(f"Image {image_path} not found.")

def gfaToFasta(gfa = "assembly.homopolymer-compressed.gfa", out_fasta=None):
    """
    Convert GFA file to FASTA format.

    Parameters
    ----------
    gfa
        Path to the GFA file. [default: assembly.homopolymer-compressed.gfa]
    out_fasta
        Path to the output FASTA file. If not provided, the output file will be named automatically. [default: None]
    
    Returns
    -------
        output FASTA file
    """
    print(f"Converting {gfa} to FASTA format")
    if out_fasta is None:
        out_fasta = gfa.replace(".gfa", ".fasta")
    if not os.path.exists(gfa):
        print(f"{gfa} not found")
        return
    if os.path.exists(out_fasta):
        print(f"{out_fasta} already exists")
        return
    with open(gfa, 'r') as f:
        with open(out_fasta, 'w') as o:
            for line in f:
                if line[0] == "S":
                    line = line.strip().split("\t")
                    o.write(f">{line[1]}\n{line[2]}\n")
    print(f"Output written to {out_fasta}")
    
def mapBetweenNodes(ref="assembly.homopolymer-compressed.fasta", query='assembly.homopolymer-compressed.fasta', threads=1, out=None, showOnly=False,
                working_directory="chromosome_assignment"):
    """
    Map the query sequences to the reference sequences using mashmap.

    Parameters
    ----------
    ref
        Path to the reference fasta file. [default: assembly.homopolymer-compressed.fasta]
    query
        Path to the query fasta file. [default: assembly.homopolymer-compressed.fasta]
    threads
        Number of threads to use. [default: 1]
    out
        Path to the output file. If not provided, the output file will be named automatically. [default: None]
    showOnly
        If set to True, the script will not be executed; it
        will only display the intended operations. [default: False]
    working_directory
        Path to the working directory. [default: chromosome_assignment]
    """
    # Ensure paths are absolute
    
    working_dir = os.path.abspath(working_directory)
    print(f"aligning {ref} to {query}")

    if not os.path.exists(working_dir):
        print(f"Creating working directory: {working_dir}")
        os.makedirs(working_dir)

    ref = os.path.abspath(ref)
    query = os.path.abspath(query)

    if out is None:
        ref_base = os.path.basename(ref)
        ref_base = re.sub(r"\.(fa|fasta|fq|fastq)(\.gz)?$", "", ref_base)

        query_base = os.path.basename(query)
        query_base = re.sub(r"\.(fa|fasta|fq|fastq)(\.gz)?$", "", query_base)

        out = f"{ref_base}_vs_{query_base}.mashmap.out"

    print(f"Reference: {ref}")
    print(f"Query: {query}")
    print(f"Output: {out}")

    # Check if the reference fasta file has an index file
    if not os.path.exists(f"{ref}.fai"):
        print("Indexing reference fasta file")
        cmd = f"samtools faidx -@ {threads} {ref}"
        run_shell(cmd, wkDir=working_dir, functionName="mapBetweenNodes_refIdx", longLog=False, showOnly=showOnly)
        
    # Check if the query fasta file has an index file
    if not os.path.exists(f"{query}.fai"):
        print("Indexing query fasta file")
        cmd = f"samtools faidx -@ {threads} {query}"
        run_shell(cmd, wkDir=working_dir, functionName="mapBetweenNodes_queryIdx", longLog=False, showOnly=showOnly)

    cmd = f"mashmap -r {ref} -q {query} -t {threads} --skipSelf --output {out}"
    if showOnly:
        # If showOnly is True, just display the command instead of executing it
        print(f"Command to be executed:\n{cmd}")
    else:
        if not os.path.exists(f"{out}"):
            print(f"Running mashmap to align {query} to {ref}")
            # Run the shell command to perform the operation
            run_shell(cmd, wkDir=working_dir, functionName="mapBetweenNodes_mashmap", longLog=False, showOnly=showOnly)