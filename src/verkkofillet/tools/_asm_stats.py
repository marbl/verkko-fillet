import sys
import shlex
import subprocess
import os

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def getT2T(obj, fasta="assembly.fasta", working_directory = "stats"):
    """
    Run the script to calculate assembly statistics, including contig presence, telomere regions, and gap locations.

    Parameters
    ----------
    obj
        An object that contains a .stats attribute, which should be a pandas DataFrame.
    fasta
        verkko assembly.
    working_directory
        output directory

    Return
    ------
        {working_directory}/assembly.gaps.bed
        {working_directory}assembly.t2t_ctgs
        {working_directory}assembly.t2t_scfs
        {working_directory}assembly.telomere.bed
    """
    script = os.path.abspath(os.path.join(script_path, "getT2T.sh"))  # Ensure absolute path
    working_dir = os.path.abspath(working_directory)  # Ensure absolute path for the working directory
    fasta = os.path.abspath(fasta)

    output_files = [
        "assembly.telomere.bed",
        "assembly.gaps.bed",
        "assembly.t2t_scfs",
        "assembly.t2t_ctgs"]

    # Check if all output files already exist
    if all(os.path.exists(os.path.join(working_dir, file)) for file in output_files):
        print("All output files already exist. Skipping getT2T.")
        return
        
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    
    cmd = f"bash {shlex.quote(script)} {shlex.quote(fasta)}"
    
    try:
        subprocess.run(
            cmd,
            stdout=subprocess.PIPE,  # Capture stdout for debugging
            stderr=subprocess.PIPE,  # Capture stderr for debugging
            shell=True,  # Allowing shell-specific syntax
            check=True,  # Raises CalledProcessError for non-zero return code
            cwd=working_dir  # Run command in the specified working directory
        )
        print("getT2T was done!")
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
        print(f"Standard output: {e.stdout.decode().strip()}")
