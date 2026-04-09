import sys
import shlex
import subprocess
import os
import pandas as pd
import re
from tqdm import tqdm


from .._run_shell import run_shell

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def detect_internal_telomere(obj, 
                             working_directory = "internal_telomere",
                             fasta = 'assembly.fasta', 
                             showOnly = False, longLog = False,
                            name= None):
    """\
    Detect internal telomere in the assembly.fasta file.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    working_directory
        The directory to store the output files. Default is "internal_telomere".
    fasta  
        The path to the assembly.fasta file. Default is "assembly.fasta".
    showOnly
        If True, the command will be printed but not executed. Default is False.
    longLog
        If True, the full output will be printed. Default is False.
    name
        The name of the output files. Default is None, which will be generated based on the input file name.

    Returns
    -------
        assembly.windows.0.5.bed
    """
    print("Starting detecting internal telomere in the assembly.fasta")

    script = os.path.abspath(os.path.join(script_path, "vgp-assembly", "telomere","telomere_analysis.sh"))
    asm = os.path.abspath(fasta)
    working_dir = os.path.abspath(working_directory)
    if name== None:
        name = name = re.sub(r"\.fasta\.gz$|\.fasta$", "", os.path.basename(asm)) + "_1"  
    
    # Define expected output files
    outputFile = os.path.join(working_dir,name, "assembly.windows.0.5.bed")

    # Check if all output files already exist
    if os.path.exists(outputFile):
        print("All output files already exist. Skipping internal telomere detection.")
        return

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    # Ensure script and working directory exist
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    cmd = f"sh {shlex.quote(script)} {name} 0.5 50000 {asm} {name}"
    run_shell(cmd, wkDir=working_dir, functionName = "detect_internal_telomere" ,longLog = longLog, showOnly = showOnly)



def runTrimming(obj, 
                trim_contig_dict,                
                original_fasta="assembly.fasta", 
                trim_bed="internal_telomere/assembly.fasta.trim.bed",
                output_fasta=None, 
                force = False,
                ):
    """\
    Trim the contigs based on the provided coordinates.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    trim_contig_dict
        Dictionary containing contig names and their corresponding trim coordinates.
    original_fasta
        The path to the original assembly.fasta file. Default is "assembly.fasta".
    trim_bed
        The path to the BED file containing trim coordinates. Default is "internal_telomere/assembly.fasta.trim.bed".
    output_fasta
        The path to the output trimmed FASTA file. Default is None, which will generate a new file name based on the original FASTA file.

    Returns
    -------
        The path to the trimmed FASTA file.

    """

    if output_fasta is None:  # Use 'is None' for comparison
        prefix = re.sub(r"\.gz|\.fasta", "", original_fasta)
        output_fasta = prefix + "_trimmed.fasta"
    
    if force and output_fasta and os.path.exists(output_fasta):
        print(f"Force is enabled. Removing existing file: {output_fasta}")
        os.remove(output_fasta)

    # Save the bed file
    rows = []
    for contig, intervals in trim_contig_dict.items():
        for start, end in intervals:
            rows.append({'contig': contig, 'from': start, 'to': end})

    trim_df = pd.DataFrame(rows)
    
    trim_df.to_csv(trim_bed, sep='\t', header=False, index=False)
    original_fasta = os.path.abspath(original_fasta)
    trim_bed = os.path.abspath(trim_bed)


            
    # Get working directory and script path
    working_dir = os.path.abspath(obj.verkko_fillet_dir)

    # Check if the final output file already exists
    if os.path.exists(output_fasta):
        print(f"{output_fasta} already exists. Exiting to avoid overwriting.")
        sys.exit(1)
    
    # Load chromosome list from FASTA index file
    fai_path = f"{original_fasta}.fai"
    
    if not os.path.exists(fai_path):
        print(f"Index file {fai_path} not found. Generating faidx index for {original_fasta}.")
        cmd = f"samtools faidx {original_fasta}"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating index file: {cmd}")
            print(f"Error message: {e}")
            sys.exit(1)
    
    # Read the fai file to get the list of chromosomes
    fai = pd.read_csv(fai_path, sep='\t', header=None, usecols=[0])
    chrList = list(fai[0])

    # Process each chromosome with a progress bar
    with tqdm(total=len(chrList), desc="Trim Chromosomes", ncols=80, colour="white") as pbar:
        for chromosome in chrList:
            # Build the command for trimming chromosomes
            if chromosome in trim_df['contig'].values:
                print(f"{chromosome} will be trimmed")
                from_value = str(trim_df.loc[trim_df['contig'] == chromosome, 'from'].values[0])
                to_value = str(trim_df.loc[trim_df['contig'] == chromosome, 'to'].values[0])
                
                cmd = f"samtools faidx {original_fasta} {chromosome}:{from_value}-{to_value} | sed -e '1d' | sed -e '1i >{chromosome}'>> {output_fasta}"
                print(f"{cmd}")
            else:
                cmd = f"samtools faidx {original_fasta} {chromosome} >> {output_fasta}"
            
            try:
                subprocess.run(cmd, shell=True, check=True)
                
            except subprocess.CalledProcessError as e:
                print(f"Error processing chromosome {chromosome}: {cmd}")
                print(f"Error message: {e}")
                # Optionally, you can log the error and continue, or raise an exception.
                continue
            
            # Update progress bar
            pbar.update(1)
    
    print(f"Trimming completed. Output file: {output_fasta}")

