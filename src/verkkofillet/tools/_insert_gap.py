import sys
import shlex
import pandas as pd
import subprocess
import os
import re


script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

def save_list_to_file(data_list, file_path="insertONTsupport.list"):
    """
    Saves a list to a file with one column. If the folder does not exist, it creates it.

    Parameters
    ----------
    data_list
        List of data to be saved.
    file_path
        Path to the file where the data will be saved.
    """
    # Extract directory from the file path and handle empty cases
    directory = os.path.dirname(file_path) if os.path.dirname(file_path) else '.'

    # Create the directory if it doesn't exist
    os.makedirs(directory, exist_ok=True)
    
    # Open the file in append mode and write the data
    with open(file_path, 'a') as f:
        for item in data_list:
            f.write(f"{item}\n")


def insertGap(gapid,
              split_reads,
              max_end_clip = 50000, 
              outputDir="missing_edge",
              graph="assembly.homopolymer-compressed.gfa"):
    """
    Find ONT support for Inserts a gap into the graph using split reads.

    Parameters
    ----------
    gapid
        Identifier for the gap.
    split_reads
        Pandas DataFrame containing reads information.
    outputDir
        Output directory for the results.
    alignGAF
        Path to alignment GAF file.
    graph
        Path to graph file.
    """

    # Ensure absolute paths
    cwd = os.getcwd()
    print(f"Current working directory: {cwd}")

    outputDir = os.path.abspath(outputDir)
    # alignGAF = os.path.abspath(alignGAF)
    graph = os.path.abspath(graph)
    
    # Check if the working directory exists
    os.makedirs(outputDir, exist_ok=True)

    try:
        # Extract Verkko script path
        script_path_proc = subprocess.run(
            ["verkko"], 
            text=True, 
            capture_output=True, 
            check=True
        )
        script_path_output = script_path_proc.stdout
        script_path_line = [line for line in script_path_output.splitlines() if "Verkko module path" in line]
        
        if not script_path_line:
            raise ValueError("Verkko module path not found in output.")
        
        verkko_path = script_path_line[0].split()[-1]
        script = os.path.abspath(os.path.join(verkko_path, "scripts", "insert_aln_gaps.py"))
        

        
    except (subprocess.CalledProcessError, ValueError) as e:
        script = os.path.join(script_path,"insert_aln_gaps.py")
        
    
    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return

    # print("Extracting reads...")

    # Ensure the column exists in split_reads
    if 'qname' not in split_reads.columns:
        print("Error: 'qname' column not found in input data.")
        return

    reads = list(set(split_reads['qname']))
    file_path = os.path.abspath(os.path.join(outputDir, f"{gapid}.missing_edge.ont_list.txt"))
    
    save_list_to_file(reads, file_path)
    print(f"The split reads for {gapid} were saved to {file_path}")

    subset_gaf = os.path.abspath(os.path.join(outputDir, f"{gapid}.missing_edge.gaf"))

    # Grep reads from GAF file
    # cmd_grep = f"grep -w -f {shlex.quote(file_path)} {shlex.quote(alignGAF)} > {shlex.quote(subset_gaf)}"

    # try:
    #     result = subprocess.run(
    #         cmd_grep, 
    #         stdout=subprocess.PIPE, 
    #         stderr=subprocess.PIPE, 
    #         shell=True, 
    #         check=True, 
    #         cwd=outputDir
    #     )
    # except subprocess.CalledProcessError as e:
    #     print(f"Command failed: {cmd_grep}")
    #     print(f"Error code: {e.returncode}")
    #     print(f"Error output: {e.stderr.decode().strip()}")
    #     return
    

    # replace "&lt;" to "<" and "&gt;" to ">" in the subset_gaf file
    split_reads.to_csv(subset_gaf, sep='\t', header=False, index=False)
    with open(subset_gaf, 'r') as file:
        data = file.read()
    data = data.replace('&lt;', '<').replace('&gt;', '>')
    with open(subset_gaf, 'w') as file:
        file.write(data)

    # Run Verkko gap insertion script
    patch_nogap = os.path.join(outputDir, f"patch.nogap.{gapid}.gaf")
    patch_gaf = os.path.join(outputDir, f"patch.{gapid}.gaf")
    patch_gfa = os.path.join(outputDir, f"patch.{gapid}.gfa")

    cmd_insert = f"python { shlex.quote(script)} {shlex.quote(graph)} {shlex.quote(subset_gaf)} 1 {max_end_clip} {shlex.quote(patch_nogap)} {shlex.quote(patch_gaf)} gapmanual y > {shlex.quote(patch_gfa)}"
    print(f"Running gap insertion command: {cmd_insert}")
    try:
        result = subprocess.run(
            cmd_insert, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True, 
            check=True, 
            cwd=outputDir
        )

        print(f"The gap filling was completed for {gapid}!")
        print(f"Patch GAF saved to: {patch_gaf}")
        
        # Display the final path contents
        try : 
            final_paths = set(pd.read_csv(patch_gaf, header=None, usecols=[5], sep='\t')[5])
            final_paths = list(final_paths)[0]


            final_paths_split = [t for t in re.split(r'(?=<)|(?=>)', final_paths) if t]
            final_paths_split_forgapFill = []
            for t in final_paths_split:
                # print(t)
                t_cl = re.sub(r'[<>]', '', t)
                if t.startswith('<') :
                    t_cl = t_cl + '-'
                elif t.startswith('>'):
                    t_cl = t_cl + '+'
                # print(t_cl)
                final_paths_split_forgapFill.append(t_cl)

            new_node_id = [item for item in final_paths_split_forgapFill if item.startswith("gap")][0]
            new_node_strand = '+' if new_node_id.endswith('+') else '-'
            new_node_id = re.sub(r'[+-]$', '', new_node_id)
            

            cmd = f"tail -3 missing_edge/patch.{gapid}.gfa > missing_edge/{gapid}.missing_edge.patching.gfa"
            subprocess.run(cmd, shell=True, check=True)
            os.remove(f"missing_edge/patch.{gapid}.gfa")

            for roi in [new_node_id] :
                fasta_path = f"missing_edge/{gapid}.missing_edge.patching.fasta"
                if os.path.exists(fasta_path):
                    cmd = f"sed -i 's/{roi}/{roi}_{gapid}/g' {fasta_path}"
                    subprocess.run(cmd, shell=True, check=True)
                
                cmd = f"sed -i 's/{roi}/{roi}_{gapid}/g' missing_edge/patch.nogap.{gapid}.gaf"
                subprocess.run(cmd, shell=True, check=True)     

                cmd = f"sed -i 's/{roi}/{roi}_{gapid}/g' missing_edge/patch.{gapid}.gaf"
                subprocess.run(cmd, shell=True, check=True)
                
                cmd = f"sed -i 's/{roi}/{roi}_{gapid}/g' missing_edge/{gapid}.missing_edge.patching.gfa"
                subprocess.run(cmd, shell=True, check=True)

            new_node_id = f"{new_node_id}_{gapid}{new_node_strand}"
            idx = final_paths_split_forgapFill.index([item for item in final_paths_split_forgapFill if item.startswith("gap")][0])
            final_paths_split_forgapFill[idx] = new_node_id
            final_paths_split_forgapFill = ','.join(final_paths_split_forgapFill)
            
            print("The new node id is:", new_node_id)
            print("The final path looks like:")
            print(final_paths_split_forgapFill)
            # return final_paths_split_forgapFill
        except Exception as e:
            print(f"Could not read final paths from {patch_gaf}: {e}")
            
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd_insert}")
        print(f"Error code: {e.returncode}")
        print(f"Error output: {e.stderr.decode().strip()}")
