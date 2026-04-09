import subprocess
import os
import pandas as pd
import pickle
import glob
import sys
import copy
import shutil
from tqdm import tqdm
from .._default_func import check_user_input, print_directory_tree,addHistory
from .._run_shell import run_shell
from datetime import datetime
import inspect

# --------------------------------------------------------------------------------
# Reading and Writing data files and AnnData objects
# --------------------------------------------------------------------------------

# Function to check user input
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

class FilletObj:
    def __init__(self):
        self.verkkoDir = None
        self.verkko_fillet_dir = None
        self.paths = None
        self.version = None
        self.species = None
        self.stats = None
        self.gaps = None
        self.gaf = None
        self.paths_freq = None
        self.qv = None
        self.history = None
        self.scfmap = None
        self.node = None
        self.edge = None

    
    def __repr__(self):
        attributes = vars(self)
        
        # Filter out None values (optional, or customize as needed)
        existing_attributes = {key: value for key, value in attributes.items() if value is not None}
        
        # Start forming the string for representation
        repr_str = f"{self.__class__.__name__}\n"
        
        # Add each attribute and its value to the string
        for attribute, value in existing_attributes.items():
            # Check if the value is a pandas DataFrame
            if isinstance(value, pd.DataFrame):
                value = ', '.join(value.columns.tolist())  # Convert the column names to a list
            
            repr_str += f"  {attribute}: {value}\n"
        
        return repr_str

def readNode(obj, graph = "assembly.homopolymer-compressed.noseq.gfa", color = "assembly.colors.csv",
             ont_cov = "8-hicPipeline/final_contigs/assembly.ont-coverage.csv",
             hifi_cov = "8-hicPipeline/final_contigs/assembly.hifi-coverage.csv"):
    
    ont_cov_df = None
    hifi_cov_df = None

    obj = copy.deepcopy(obj)
    if not os.path.exists(graph):
        raise FileNotFoundError(f"File {graph} not found")
    if not os.path.exists(color):
        raise FileNotFoundError(f"File {color} not found")
    
    if not os.path.exists(ont_cov):
        print(f"File {ont_cov} not found")
    else:
        ont_cov_df = pd.read_csv(ont_cov, sep='\t', header=0, usecols=[0,1])
        ont_cov_df.columns = ['node', 'ont_cov']
        ont_cov_df['ont_cov'] = pd.to_numeric(ont_cov_df['ont_cov'], errors='coerce')
    
    if not os.path.exists(hifi_cov):
        print(f"File {hifi_cov} not found")
    else:
        hifi_cov_df = pd.read_csv(hifi_cov, sep='\t', header=0, usecols=[0,1])
        hifi_cov_df.columns = ['node', 'hifi_cov']
        hifi_cov_df['hifi_cov'] = pd.to_numeric(hifi_cov_df['hifi_cov'], errors='coerce')

    print(f"Reading {graph}")
    nodeLen = pd.read_csv(graph, sep='\t', header=None)
    nodeLen = nodeLen[nodeLen[0] == "S"]
    nodeLen = nodeLen.iloc[:,[1,3]]
    nodeLen.columns = ['node', 'len']
    nodeLen['len'] = nodeLen['len'].str.replace(r'^LN:i:', '', regex=True)
    nodeLen['len'] = pd.to_numeric(nodeLen['len'], errors='coerce')  # Handle non-numeric values gracefully
    print(f"Reading {color}")
    color = pd.read_csv(color, sep='\t', header=0)
    df = pd.merge(nodeLen, color, on='node', how='left')

    if ont_cov_df is not None:
        df = pd.merge(df, ont_cov_df, on='node', how='left')
    if hifi_cov_df is not None:
        df = pd.merge(df, hifi_cov_df, on='node', how='left')

    obj.node = df
    obj = addHistory(obj, f"node file is loaded from {graph}", {inspect.currentframe().f_code.co_name})
    print(f"Node information is stored in obj.node")
    return obj

def readEdge(obj, graph = "assembly.homopolymer-compressed.noseq.gfa"):
    obj = copy.deepcopy(obj)
    nodeLen = pd.read_csv(graph, sep='\t', header=None)
    nodeLen = nodeLen[nodeLen[0] == "L"]
    nodeLen = nodeLen.loc[:,1:5]
    nodeLen.columns = ["node1", "node1_strand", "node2", "node2_strand", "overlap"]
    obj.edge = nodeLen
    obj = addHistory(obj, f"edge file is loaded from {graph}", {inspect.currentframe().f_code.co_name})
    return obj

def readPath(obj, paths_path = "assembly.paths.tsv"):
    """
    """
    if not os.path.exists(paths_path):
        raise FileNotFoundError(f"File {paths_path}")
    
    obj = copy.deepcopy(obj)

    # Load paths file if it exist
    obj.paths = pd.read_csv(paths_path, header=0, sep='\t', index_col=None)
    print("Path file loaded successfully.")
    obj = addHistory(obj, f"path file is loaded from {paths_path}", {inspect.currentframe().f_code.co_name})
    
    return obj

def readScfmap(obj, scfmap_path = "assembly.scfmap"):
    """
    """
    obj = copy.deepcopy(obj)
    if not os.path.exists(scfmap_path):
        raise FileNotFoundError(f"File {scfmap_path} not found")

    scfmap = pd.read_csv(scfmap_path, sep = ' ', header = None)
    scfmap.columns = ['info','contig','pathName']
    scfmap= scfmap.loc[scfmap['info']=='path']
    del scfmap['info']
    scfmap.reset_index(drop=True, inplace=True)
    obj.scfmap = scfmap
    print("scfmap file loaded successfully.")
    
    
    obj = addHistory(obj, f"scfmap file is loaded from {scfmap_path}", {inspect.currentframe().f_code.co_name})
    return obj



def read_Verkko(verkkoDir, 
                verkko_fillet_dir=None, 
                paths_path="assembly.paths.tsv",
                force = False,
                scfmap_path = "assembly.scfmap", 
                version=None, 
                species=None, 
                graph = "assembly.homopolymer-compressed.noseq.gfa",
                color = "assembly.colors.csv",
                lock_original_folder = True, showOnly = False, longLog = False):
    """
    Prepares the Verkko environment by creating necessary directories, locking the original directory, 
    and loading the paths file for further processing.

    Parameters
    ----------
    verkkoDir
        Base directory of Verkko data.
    verkko_fillet_dir
        Target directory for fillet data. Defaults to None.
    paths_path 
        Path to 'assembly.paths.tsv' file. Defaults to 'assembly.paths.tsv'.
    version
        Version of the data. Defaults to None.
    species
        Species name. Defaults to None.
    lock_original_folder
        Whether to lock the original directory. Defaults to True.

    Returns
    -------
    obj : FilletObj
        A FilletObj instance with the configured directories and loaded paths data.
    """
    # make filletObj
    obj = FilletObj()
    
    verkkoDir = os.path.realpath(verkkoDir)
    
    # set verkko_fillet output dir
    if verkko_fillet_dir == None:
        verkko_fillet_dir = os.path.join(verkkoDir+"_verkko_fillet")

    # check the verkko fillet output dir
    if os.path.exists(verkko_fillet_dir) and force == False:
        print(f"The Verkko fillet target directory already exists: {verkko_fillet_dir}")
        print(f"If you didn't mean this, please set another directory or for overwirting, please use force= True")
    else:
        # Create the directory if it does not exist
        print(f"The Verkko fillet target directory has been created and set to: {verkko_fillet_dir}")
        print("All temporary and output files will be written to this directory.")
        script = os.path.abspath(os.path.join(script_path, "make_verkko_fillet_dir.sh"))
        cmd=f"sh {script} {verkkoDir} {verkko_fillet_dir}"
        run_shell(cmd, wkDir=verkko_fillet_dir, functionName = "make_verkko_fillet_dir" ,longLog = longLog, showOnly = showOnly)
        
    
    # lock original verkko folder to prevent mess up
    if lock_original_folder :
        print(f"Lock the original Verkko folder to prevent it from being modified.")
        script = os.path.abspath(os.path.join(script_path, "lock_folder.sh"))
        cmd=f"sh {script} {verkkoDir}"
        run_shell(cmd, wkDir=verkko_fillet_dir, functionName = "lock_original_folder" ,longLog = longLog, showOnly = showOnly)

    # Set the additional attributes on the object
    obj.species = species
    obj.verkkoDir = verkkoDir
    obj.verkko_fillet_dir = verkko_fillet_dir
    obj.version = version
    obj.history = pd.DataFrame({"timestamp": [datetime.now()],"activity": [f"verkko-fillet obj is generated. from : {verkkoDir}, outdir : {verkko_fillet_dir}"], "function" : "read_Verkko"})
    # obj = addHistory(obj, f"verkko-fillet obj is generated. from : {verkkoDir}", {inspect.currentframe().f_code.co_name})
    
    print(f"change working directory: {verkko_fillet_dir}")
    os.chdir(verkko_fillet_dir)
    
    # read Path file
    if paths_path != None:
        obj = readPath(obj, paths_path)

    # read scfmap file
    if scfmap_path != None:
        obj = readScfmap(obj, scfmap_path)

    # read node file 
    if graph != None:
        obj = readNode(obj, graph, color)
        obj = readEdge(obj, graph)
    
    return obj

def save_Verkko(obj,
                fileName: str):
    """\
    Save the Verkko fillet object to a file using pickle.

    Parameters
    ----------
    obj
        The Verkko fillet object to be saved.
    fileName
        The name of the file to save the object to.
    """
    print("save verkko fllet obj to -> " + fileName)
    obj = addHistory(obj,f"Writing verkko-fillet obj to {fileName}", 'save_Verkko')
    with open(fileName, "wb") as f:
        pickle.dump(obj, f)

def load_Verkko(fileName):
    """\
    Load the Verkko fillet object from a file using pickle.

    Parameters
    ----------
    fileName
        The name of the file to load the object from.

    Returns
    -------
    obj: object
        The loaded Verkko fillet object.
    """
    print("load verkko fllet obj from <- " + fileName)
    # Open the file in read-binary mode
    with open(fileName, "rb") as f:
        # Load the object from the file using pickle
        obj = pickle.load(f)
        
    obj = addHistory(obj,f"Reading verkko-fillet obj from {fileName}", 'load_Verkko')
    
    return obj
    
def hard_copy_symlink(symlink_path, destination_path):
    """
    Creates a hard copy of the file pointed to by the symbolic link.
    
    Parameters
    ----------
    symlink_path : str
        The path to the symbolic link.
    destination_path : str
        The path to the destination where the hard copy will be created.
    """
    if os.path.islink(symlink_path):
        # Get the target of the symbolic link
        target_path = os.readlink(symlink_path)
        #print(f"Symbolic link points to: {target_path}")

        # Copy the actual file to the destination
        shutil.copy(target_path, destination_path)
        #print(f"Hard copy of the symlink created at: {destination_path}")
    else :
        shutil.copy(symlink_path, destination_path)



def updateCNSdir_missingEdges(obj, new_folder_path, 
                              missing_edge_dir ="missing_edge",
                              final_gaf = "assembly.fixed.paths.gaf",  showOnly = False, longLog = False):
    """
    Updates the CNS directory by handling missing edges and creating necessary symbolic links or files.
    
    Parameters
    ----------
    obj
        Object containing the original verkko directory path.
    new_folder_path
        Path to the new folder to be updated.

    Returns
    -------
        new folder with updated files and symbolic links for missing edges
    """
    newFolder = os.path.abspath(new_folder_path)
    filletDir = os.path.abspath(obj.verkko_fillet_dir)  # Define oriDir only once
    verkkoDir = os.path.abspath(obj.verkkoDir)
    missing_edge_dir = os.path.abspath(missing_edge_dir)
    
    # Check if the new folder exists
    if not os.path.exists(newFolder):
        print("New verkko folder for CNS is not exists!")
        return
    if not os.path.exists(missing_edge_dir):
        print("New verkko folder for CNS is not exists!")
        return
    
    file_list = glob.glob(os.path.join(missing_edge_dir, "patch.gapid_*.gaf"))  # Check if there are files in the missing edge directory
    # get gapid_* from the file name and check if there are corresponding .gfa files
    gapid_list = [os.path.basename(file).split(".")[1] for file in file_list]
    print(f"Found gapid list that have been filled by missing edges function : {gapid_list}")
    script = os.path.abspath(os.path.join(script_path, "_updateCNSdir_missingEdges.sh"))
    cmd=f"sh {script} {filletDir} {verkkoDir} {newFolder} {final_gaf} {missing_edge_dir}"
    run_shell(cmd, wkDir=filletDir, functionName = "make_verkko_fillet_dir" ,longLog = longLog, showOnly = showOnly)

def checkFiles(folder):
    files = [
        "1-buildGraph", "2-processGraph", "3-align", "3-alignTips", "4-processONT", "5-untip", 
        '6-layoutContigs/combined-alignments.gaf', '6-layoutContigs/combined-edges.gfa',
        '6-layoutContigs/consensus_paths.txt', '6-layoutContigs/nodelens.txt',
        '6-layoutContigs/ont-gapfill.txt', '6-layoutContigs/ont.alignments.gaf',  # Removed trailing space
        '7-consensus/ont_subset.fasta.gz', '7-consensus/ont_subset.id',
        'hifi-corrected.fasta.gz', "6-layoutContigs/combined-nodemap.txt", "6-layoutContigs/hifi.alignments.gaf",
    ]
    notExist = []
    for file in tqdm(files):
        if not os.path.exists(os.path.join(folder, file)):  # Fixed logic
            notExist.append(file)
    
    if notExist:
        print(f"Following files are missing in the folder: {notExist}")
    else:
        print("All files exist.")



def mkCNSdir(obj, new_folder_path, final_gaf = "assembly.fixed.paths.gaf", 
             missingEdge = False, missing_edge_dir = "missing_edge",):
    """\
    Creates a new CNS directory by creating symbolic links to the original verkko directory.

    Parameters
    ----------
    obj
        Object containing the original verkko directory path.
    new_folder_path
        Path to the new folder to be created.
    final_gaf
        Path to the final GAF file. Default is "final_rukki_fixed.paths.gaf".
    missingEdge
        Whether to handle missing edges. Default is False.
    tmp_id
        Path to the temporary ID file. Default is "missing_edge/ont_subset.tmp.id".
    tmp_fasta
        Path to the temporary FASTA file. Default is "missing_edge/ont_subset.tmp.fasta".

    Returns
    -------
        new folder with mendatory files and symbolic links
    """
    print(f"Creating a new verkko folder for CNS at: {new_folder_path}")
    print(f"Copying the final GAF file from: {final_gaf}")

    # if missing_edge_dir is exist, set missingEdge to True
    if os.path.exists(missing_edge_dir):
        missingEdge = True
        print(f"missing_edge_dir {missing_edge_dir} is exist, set missingEdge to True")
        print(f"missingEdge mode is set to: {missingEdge}")
        
    newFolder = os.path.abspath(new_folder_path)
    verkko_fillet_dir = os.path.abspath(obj.verkko_fillet_dir)  # Define original directory
    verkkoDir = os.path.abspath(obj.verkkoDir)  # Define original directory
    final_gaf = os.path.abspath(final_gaf)  # Define final GAF file
    missing_edge_dir = os.path.abspath(os.path.join(verkko_fillet_dir, missing_edge_dir))
    # Create the new folder (only if it doesn't exist)
    os.makedirs(newFolder, exist_ok=True)

    # Check if the folder already exists
    if os.path.exists(newFolder) and os.listdir(newFolder):
        print("New verkko folder for CNS already exists and is not empty!")
        return

    # Create symbolic links for directories
    for folder in ["1-buildGraph", "2-processGraph", "3-align", "3-alignTips", "4-processONT", "6-rukki", "5-untip", 
                   "hifi-corrected.fasta.gz", "assembly.colors.csv", "assembly.homopolymer-compressed.noseq.gfa",]:
        source_path = os.path.join(verkkoDir, folder)
        link_path = os.path.join(newFolder, folder)

        if not os.path.exists(source_path):
            print(f"Warning: Directory or file{folder} does not exist in the original Verkko directory.")
            continue

        # Check if symlink already exists
        if not os.path.exists(link_path):
            os.symlink(source_path, link_path)

    if missingEdge:
        print(f"Handling missing edges in the new folder: {newFolder}")
        updateCNSdir_missingEdges(obj = obj,  new_folder_path= new_folder_path, final_gaf = final_gaf, missing_edge_dir =missing_edge_dir)
        
    else:
        # Create additional folders
        os.makedirs(os.path.join(newFolder, "6-layoutContigs"), exist_ok=True)
        os.makedirs(os.path.join(newFolder, "7-consensus"), exist_ok=True)

        # Create symbolic links for specific files
        files_to_link = [
            "6-layoutContigs/combined-nodemap.txt",
            "6-layoutContigs/hifi.alignments.gaf",
            "6-layoutContigs/combined-alignments.gaf",
            "6-layoutContigs/combined-edges.gfa",
            "6-layoutContigs/nodelens.txt",
            "7-consensus/ont_subset.fasta.gz",
            "7-consensus/ont_subset.id",
            "6-layoutContigs/ont-gapfill.txt",
            "6-layoutContigs/ont.alignments.gaf",
        ]

        for file in files_to_link:
            source_file = os.path.join(verkkoDir, file)
            link_file = os.path.join(newFolder, file)

            if not os.path.exists(source_file):
                print(f"Warning: File {file} does not exist in the original Verkko directory.")
                continue

            if not os.path.exists(link_file):
                os.symlink(source_file, link_file)

    # Ensure `final_gaf` exists before copying
    final_gaf = os.path.abspath(final_gaf)
    target_gaf = os.path.join(newFolder, "6-layoutContigs", "consensus_paths.txt")
    

    if not os.path.exists(final_gaf):
        print(f"Warning: File {final_gaf} does not exist in the original Verkko directory.")
    else:
        subprocess.run(["cp", final_gaf, target_gaf], check=True)

    print("✅ All files are updated! The new folder is ready for verkko-cns.")
    print(" ")
    print(f"Checking the new folder for missing files...")
    checkFiles(newFolder)


testDir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../data/'))
def loadGiraffe():
    """\
    Load the object of Giraffe genome from a file using pickle.

    Returns
    -------
    obj : object
        The loaded Giraffe genome object.
    """
    fileName = f"{testDir}/test_giraffe/giraffe_before_gap_filling.pkl"
    obj = load_Verkko(fileName)
    return obj