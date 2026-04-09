import pandas as pd
import logging
import re
import copy
import time
from tqdm import tqdm
from .._default_func import addHistory
# Configure logging
logging.basicConfig(level=logging.INFO)

def path_to_gaf(input_string):
    """
    Converts a path string into a GAF format.

    Parameters
    ----------
        input_string (str): The input path string to be converted.

    Returns
    -------
        str: The converted GAF format string.
    """
    
    # Split the input string by commas
    items = input_string.split(',')
    
    # Process each item
    result = []
    for item in items:
        if item.endswith('+'):
            result.append('>' + item[:-1])
        elif item.endswith('-'):
            result.append('<' + item[:-1])
        elif item.endswith(']'):
            result.append(item)
        else :
            error = "Invalid input string: " + input_string
            print(error)
            return

    
    # Join the processed items into a single string
    return ''.join(result)


def progress_bar(current, total):
    """
    Displays a progress bar in the console.
    
    Parameters
    ----------
        current (int): Current progress.
        total (int): Total progress.
    """
    progress = int((current / total) * 50)
    bar = "[" + "=" * progress + " " * (50 - progress) + "]"
    print(f"\r{bar} {current}/{total} gaps filled", end="")
    print(" ")

def checkGapFilling(obj):
    """
    This function checks and prints the number of filled gaps in the 'gap' DataFrame
    and shows the progress bar for gap filling.

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'gap' DataFrame in obj.gaps.
    
    Returns
    ----------
        The updated 'gap' DataFrame.
    """
    obj_sub = copy.deepcopy(obj)
    gap = obj_sub.gaps  # Assuming gap is the DataFrame containing gap information
    
    gap = gap.loc[~gap['cat'].isin(['addContig','connectContig'])]
    total = gap.shape[0]  # Total number of gaps
    # gap['finalGaf'] = gap['finalGaf'].str.replace('<', '&lt;').str.replace('>', '&gt;')
    # gap['done'] = gap['fixedPath'].apply(lambda x: "✅" if x else "")
    # Count the number of non-empty 'fixedPath' entries
    current = gap['fixedPath'].apply(lambda x: pd.notna(x) and x != "").sum()
    
    # Print the current and total number of filled gaps
    # print(f"Number of filled gaps: {current} of total gaps: {total}")

    # Call the progress_bar function to show the filling progress
    progress_bar(current, total)
    
    return gap

def transform_path(elements):
    """
    Transforms elements of the path for gap filling.

    Parameters
    ----------
    elements
        A list of elements in the path.
    
    Returns
    -------
    list
        A list of transformed elements.
    """
    return [
        (">" + elem[:-1] if elem.endswith("+") else "<" + elem[:-1]) if not elem.startswith("[") else elem
        for elem in elements
    ]

def check_match(gap_value, element, position):
    """
    Checks if a specific gap matches the given element.
    
    Parameters
    ----------
    gap_value
        The gap value from the DataFrame.
    element
        The element to match.
    position
        The position in the gap (0 for start, 2 for end).
    
    Returns
    -------
        "match" if matches, else "notMatch".
    """
    return "match" if gap_value[position] == element else "notMatch"

def fillGaps(obj, gapId, final_path, notes = None, cat = "gapFill_with_evidence"):
    """
    Fills gaps for a specific gapId, updates the 'fixedPath', 'startMatch', 'endMatch', and 'finalGaf' columns.
    
    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'gap' DataFrame in obj.gaps.
    gapId
        The identifier for the gap.
    final_path
        The final path to fill the gap.
    notes
        The notes to add to the gap. Default is None.
    cat 
        "gapFill_with_evidence", "random_assign", "maynot_correct"


    Returns
    -------
    obj : object
        The updated verkko fillet object.
    """
    obj_sub = copy.deepcopy(obj)
    gap = obj_sub.gaps  # The DataFrame containing gap data

    # Ensure the gapId exists
    if gapId not in gap['gapId'].values:
        raise ValueError(f"gapId {gapId} not found in the DataFrame.")

    # Handle empty final_path
    if final_path == "":
        gap.loc[gap['gapId'] == gapId, ['fixedPath', 'startMatch', 'endMatch', 'finalGaf']] = ""
        print(f"gapId {gapId}: 'final_path' is empty. Other columns have been reset to 'NA'.")
    else:
        # Update the 'fixedPath' column for the matching gapId
        gap.loc[gap['gapId'] == gapId, 'fixedPath'] = final_path

        elements = final_path.replace(" ", "").split(",")
        modified_elements = transform_path(elements)
        modified_path = "".join(modified_elements)
        print(f"final path : {modified_path}")

        # Update the 'finalGaf' column for the matching gapId
        gap.loc[gap['gapId'] == gapId, 'finalGaf'] = modified_path

        # Retrieve the matching row for further updates
        gap_row = gap.loc[gap['gapId'] == gapId].iloc[0]

        # Check the direction and update 'startMatch' and 'endMatch'
        gap.loc[gap['gapId'] == gapId, 'startMatch'] = check_match(gap_row.gaps, elements[0], 0)
        gap.loc[gap['gapId'] == gapId, 'endMatch'] = check_match(gap_row.gaps, elements[-1], -1)
        
        if notes is not None:
            gap.loc[gap['gapId'] == gapId, 'notes'] = notes
        if cat is not None:
            if cat not in ["gapFill_with_evidence", "random_assign", "maynot_correct"]:
                raise ValueError(f"cat {cat} not found in the DataFrame.")
            gap.loc[gap['gapId'] == gapId, 'cat'] = cat

        print(f"Updated gapId {gapId}!")
        print(" ")
        if check_match(gap_row.gaps, elements[0], 0) == "match" :
            print("✅ The start node and its direction match the original node.")
        else :
            print("❌ The start node and its direction do not match the original node.")
        
        if check_match(gap_row.gaps, elements[-1], -1) == "match" :
            print("✅ The end node and its direction match the original node.")
        else :
            print("❌ The end node and its direction do not match the original node.")
        
    # Count remaining empty strings or 'NA' in 'finalGaf
    obj_sub.gaps = gap
    
    obj_sub = addHistory(obj_sub, f"{gapId} filled with {final_path}", 'fillGaps')
    # Show progress after each gap filled
    checkGapFilling(obj_sub)
    
    # Return the updated object
    
    return obj_sub

# Reset the index of the 'gap' DataFrame


def preprocess_path(path_str):
    path_str = path_str.replace("\[", ",\[")
    path_str = path_str.replace("\]", "\],")
    split_path = re.split(r',', path_str)
    return [p for p in split_path if p.strip()]  # Remove empty elements

def connectContigs(obj, contig, contig_to,  at = "left", gap = "[N5000N:connectContig]", flip = False):
    """
    Connects two contigs by adding a gap between them.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used. 
    contig
        The name of the contig to connect. 
    contig_to
        The name of the contig to connect to. 
    at
        The position to connect the contig. Default is "left". Other option is "right".
    gap
        The gap to add between the contigs. Default is "[N5000N:connectContig]". 
    flip
        Whether to flip the path from contig. Default is False. If True, the path will be flipped.

    Returns
    -------
        The updated VerkkoFillet object with new gap added.
    """
    obj_sub = copy.deepcopy(obj)
    pathdb = obj_sub.paths.copy()
    gapdb = obj_sub.gaps.copy()
    gap= gap.replace(" ", "")
    path1_raw = pathdb.loc[pathdb['name'] == contig_to]["path"].values
    path2_raw = pathdb.loc[pathdb['name'] == contig]["path"].values

    path2_path = preprocess_path(path2_raw[0])

    gapid_add= f"gapid_{str(gapdb['gapId'].str.replace('gapid_','').astype(int).max()+1)}"

    if flip:
        path2_path = path2_path[::-1]
        path2_path = [s.translate(str.maketrans("+-", "-+")) for s in path2_path]

    if at == "left":
        marker = ["startMarker"]
        fixedPath = path2_path + [gap]+ marker
    if at == "right":
        marker = ['endMarker']
        fixedPath = marker + [gap] +  path2_path
        
    fixedPath = ','.join(fixedPath)
    gap_new_line = pd.DataFrame({"gapId": gapid_add, 
                  "name" : [contig],
                  "gaps" : marker,
                  "notes" : f"connected {contig} to {contig_to} at {at} with flip {flip} {gap}",
                  "fixedPath": fixedPath,
                  "startMatch" : "",
                  "endMatch" : "",
                  "finalGaf" : "",
                  "done" : True,
                  "cat" : "connectContig",
                  })

    gapdb = pd.concat([gapdb,gap_new_line], ignore_index=True)
    obj_sub = addHistory(obj_sub, f"{gapid_add} was created", 'connectContig')
    
    obj_sub.gaps = gapdb

    print(f"Connected {contig} to {contig_to} at {at} with flip {flip}")
    print(f"{contig} was merged to {contig_to} in obj.paths")
    print(f"{contig} was replaced with {contig_to} in obj.gaps")
    print(f"New gap was added to obj.gaps with gapId {gapid_add}")
    return obj_sub

def addContigs(obj, new_path, contig_to,  at = "left", gap = "[N5000N:connectContig]"):
    """
    Connects two contigs by adding a gap between them.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used. 
    new_path
        The new path to connect. 
    contig_to
        The name of the contig to connect to. 
    at
        The position to connect the contig. Default is "left". Other option is "right".
    gap
        The gap to add between the contigs. Default is "[N5000N:connectContig]". 
    flip
        Whether to flip the path from contig. Default is False. If True, the path will be flipped.

    Returns
    -------
        The updated VerkkoFillet object with new gap added.
    """
    obj_sub = copy.deepcopy(obj)
    pathdb = obj_sub.paths.copy()
    gapdb = obj_sub.gaps.copy()
    gap= gap.replace(" ", "")

    path2_path = new_path

    gapid_add= f"gapid_{str(gapdb['gapId'].str.replace('gapid_','').astype(int).max()+1)}"

    if at == "left":
        marker = ["startMarker"]
        if gap != "":
            fixedPath = [path2_path] + [gap]
        else:
            fixedPath = [path2_path]
    if at == "right":
        marker = ['endMarker']
        if gap != "":
            fixedPath =[gap] +  [path2_path]
        else:
            fixedPath = [path2_path]
        
    fixedPath = ','.join(fixedPath)
    gap_new_line = pd.DataFrame({"gapId": gapid_add, 
                  "name" : [contig_to],
                  "gaps" : marker,
                  "notes" : f"connected new path to {contig_to} at {at} with {gap}",
                  "fixedPath": fixedPath,
                  "startMatch" : "",
                  "endMatch" : "",
                  "finalGaf" : "",
                  "done" : True,
                  "cat" : "addContig",})

    gapdb = pd.concat([gapdb,gap_new_line], ignore_index=True)
    obj_sub = addHistory(obj_sub, f"{gapid_add} was created", 'addContig')
    
    obj_sub.gaps = gapdb

    return obj_sub

def deleteGap(obj, gapId):

    """
    Deletes a gap from the 'gap' DataFrame for a specific gapId.

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'gap' DataFrame in obj.gaps.
    gapId
        The identifier for the gap to delete.

    Returns
    -------
        The updated verkko fillet object.
    """
    obj_sub = copy.deepcopy(obj)
    if gapId not in obj.gaps['gapId'].values:
        raise ValueError(f"gapId {gapId} not found in the DataFrame.")
    gaps = obj_sub.gaps.copy()
    gaps = gaps.loc[gaps['gapId'] != gapId]
    obj_sub.gaps = gaps
    obj_sub = addHistory(obj_sub, f"{gapId} was removed", 'deleteGap')
    return obj_sub

def saveGapNodes(obj, save = "gapNodes.tsv"):

    gapdf = obj.gaps.copy()

    gapdf=gapdf.loc[(gapdf['cat'] != "") & (pd.notna(gapdf['cat']))]
    gapdf = gapdf[['cat', 'fixedPath']]

    gapdf['fixedPath'] = gapdf['fixedPath'].str.split(',\s*')  # split on comma + optional space

    # Step 2: Explode to create new rows for each element
    gapdf = gapdf.explode('fixedPath')
    gapdf['fixedPath'] = gapdf['fixedPath'].replace(r'-$', '', regex=True)
    gapdf['fixedPath'] = gapdf['fixedPath'].replace(r'\+$', '', regex=True)
    gapdf = gapdf[~gapdf['fixedPath'].str.contains(r'\[', na=False)]

    gapdf['color'] = ""
    gapdf.loc[gapdf['cat'] == 'gapFill_with_evidence', 'color'] = '#4ceb34' #green
    gapdf.loc[gapdf['cat'] == 'maynot_correct', 'color'] = '#eb34e5' #red
    gapdf.loc[gapdf['cat'] == 'random_assign', 'color'] = '#a020f0' #purple
    gapdf.loc[gapdf['cat'] == 'connectContig', 'color'] = '#fff34f' #yellow
    gapdf.loc[gapdf['cat'] == 'addContig', 'color'] = '#fff34f' #yellow
    gapdf.columns = ['cat', 'node', 'color']
    gapdf = gapdf[['node','cat','color']]
    print("Saving gap nodes to ", save)
    gapdf.to_csv(save, index=False, sep='\t')

def writeFixedPaths(obj, save_path = "assembly.fixed.paths.tsv", save_gaf = "assembly.fixed.paths.gaf", save_gapNode= "gap.nodes.tsv"):
    """
    Writes the fixed paths to a file.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used. 
    save_path
        The file path to save the fixed paths. Default is "assembly.fixed.paths.tsv".
    save_gaf
        The file path to save the fixed gaf. Default is "assembly.fixed.paths.gaf".

    Returns
    -------
        The fixed paths and gaf saved to the specified file paths.
    """
    print(f"Starting to write fixed paths to {save_path} and {save_gaf}...")
    print(" ")
    obj = copy.deepcopy(obj)
    path = obj.paths.copy()
    gaps = obj.gaps.copy()
    gaps.reset_index(drop=True, inplace=True)
    path.reset_index(drop=True, inplace=True)
    gaps["fixedPath"] = gaps["fixedPath"].fillna("empty")
    gaps["fixedPath"] = gaps["fixedPath"].replace("", "empty")
    gaps["fixedPath"] = gaps["fixedPath"].replace(" ", "empty")
    # show stats
    statDb = path.groupby('rm').size().reset_index(name='counts').sort_values(by='counts', ascending=False)
    cat = list(statDb['rm'])
    count = list(statDb['counts']) 
    print(f"The total number of original paths is {len(path)}")
    print(f"The number of paths that were removed:")
    for i in range(len(cat)):
        print(f"{cat[i]}: {count[i]}")
    print(" ")

    # Fix the paths with gap filling
    for i in tqdm(range(len(gaps)), desc="Fixing paths"):
        fixed_gap = gaps.loc[i, "fixedPath"].replace(" ", "")
        ori_gap = gaps.loc[i, "gaps"]
        name = gaps.loc[i, "name"]
        gapId = gaps.loc[i, "gapId"]
        if ori_gap not in ["startMarker", "endMarker"] and fixed_gap != "empty":
            if isinstance(ori_gap, list):
                ori_gap = ','.join(ori_gap)
        
            #if name == "sire_compressed.k31.hapmer_from_utig4-1859":
                #print(f"gapid_31: {ori_gap} {fixed_gap}")

            ori_path = path.loc[path['name'] == name, "path"].values[0].replace(" ", "")
            
            if ori_gap in ori_path:
                ori_path = ori_path.replace(ori_gap, fixed_gap)
                
            else:
                print(f"{gapId} {name} {ori_gap} not in original path")
                print(f"{gapId} trying to fix with new start and end nodes in the final path")
                ori_gap[0] = fixed_gap.replace(" ", "").split(",")[0] # split the path into a list
                ori_gap[-1] = fixed_gap.replace(" ", "").split(",")[-1] # split the path into a list
                
            
            # print(fixed_path)
            fixed_path = ori_path.replace(ori_gap, fixed_gap).replace(" ", "")
            path.loc[path['name'] == name, "path"] = fixed_path
    
    # Fix the connecting paths
    # gaps_connect = gaps[gaps['gaps'].isin(['startMarker', 'endMarker'])].reset_index(drop=True)
    gaps_connect = gaps[gaps['cat'].isin(['connectContig'])].reset_index(drop=True)

    for i in range(len(gaps_connect)):
        note= gaps_connect.loc[i, 'notes'].split(" ")
        marker = gaps_connect.loc[i, 'gaps']
        
        target = note[3]
        source = note[1]
        flip = note[8]
        gap = note[9]
        
        print(f"{marker} {target} {source} {flip}")

        target_path = path.loc[path['name'] == target, 'path'].values[0]
        source_path = path.loc[path['name'] == source, 'path'].values[0]

        if flip == "True":
            source_path = source_path[::-1]
            source_path = [s.translate(str.maketrans("+-", "-+")) for s in source_path]
        
        source_path = source_path.replace(" ", "")
        target_path = target_path.replace(" ", "")
        
        if marker == "startMarker":
            newPath = [source_path] + [gap] + [target_path]
            newPath = [item for item in newPath if item != ""]
            newPath = ','.join(newPath)
            newPath = newPath.replace("startMarker", "")
            
        if marker == "endMarker":
            newPath = [target_path] + [gap] + [source_path]
            newPath = [item for item in newPath if item != ""]
            newPath = ','.join(newPath)
            newPath = newPath.replace("endMarker", "")
            
        newPath = newPath.replace(" ", "")
        path.loc[path['name'] == target, 'path'] = newPath


    gaps_connect = gaps[gaps['cat'].isin(['addContig'])].reset_index(drop=True)

    for i in range(len(gaps_connect)):
        marker = gaps_connect.loc[i, 'gaps']
        fixed_path = gaps_connect.loc[i, 'fixedPath'].replace(" ", "")
        target = gaps_connect.loc[i, 'name']
        target_path = path.loc[path['name'] == target, 'path'].values[0].replace(" ", "")

        if marker == "startMarker":
            newPath = [fixed_path] + [target_path]
            newPath = [item for item in newPath if item != ""]
            newPath = ','.join(newPath)
            newPath = newPath.replace("startMarker", "")

        elif marker == "endMarker":
            newPath = [target_path] + [fixed_path]
            newPath = [item for item in newPath if item != ""]
            newPath = ','.join(newPath)
            newPath = newPath.replace("endMarker", "")

        newPath = newPath.replace(" ", "")
        path.loc[path['name'] == target, 'path'] = newPath

    print(f"fix the connecting path is done!")
    
    # Remove paths that are covered by others
    excludelst = []

    for i in tqdm(range(len(path)), ncols=80, colour="white", desc="Excluding duplicated paths"):
        path_list = path.loc[i, 'path'].split(",")
        path_list = [re.sub(r"[+-]$", "", s) for s in path_list]
        path_name = path.loc[i, 'name']

        # Get all other paths (excluding current row), then flatten into one list
        path_etc = path.drop(i)['path'].str.split(",").explode().tolist()
        path_etc = [re.sub(r"[+-]$", "", s) for s in path_etc]

#        if all(item in path_etc for item in path_list):
#            excludelst.append(path_name)
        if set(path_list).issubset(set(path_etc)):
            excludelst.append(path_name)


    path.loc[path['name'].isin(excludelst), 'rm'] = "rm_covered_by_others"

    pat_filtered = path.loc[path['rm'].str.startswith('keep')]

    
    if "index" in pat_filtered.columns:
        pat_filtered = path.drop(columns = "index")
    #if "rm" in path.columns:
    #    path = path.drop(columns = "rm")
    print(f"The number of paths that were kept: {len(pat_filtered)}")
    pat_filtered = pat_filtered.reset_index(drop = True)

    gaf = pat_filtered.copy()
    gaf['path'] = pat_filtered['path'].apply(path_to_gaf).apply(lambda x: x.replace(' ', '') if x is not None else '')
    print(" ")
    print(f"The total number of final paths is {len(pat_filtered)}")
    
    pat_filtered.to_csv(save_path, sep = "\t", index = False)
    print(f"Fixed paths were saved to {save_path}")
    
    gaf.to_csv(save_gaf, sep = "\t", index = False)
    print(f"Fixed gaf were saved to {save_gaf}")

    path.loc[path['path'].str.contains(r'\[', regex=True),'name'] 
    contigs_with_gaps = path.loc[path['path'].str.contains(r'\[', regex=True), 'name'].tolist()
    print(" ")
    print("the contigs that have gaps:", contigs_with_gaps)


    saveGapNodes(obj, save = save_gapNode)
    return path

def checkDisconnectNode(obj, min_hpc_len = 100_000):
    """
    This function checks for disconnected nodes in the graph and filters them based on a minimum length.

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'node' and 'edge' DataFrames in obj.node and obj.edge.
    min_hpc_len
        The minimum length for filtering out disconnected nodes. Default is 100,000.

    Returns
    -------
        The updated verkko fillet object with disconnected nodes filtered out.
    """
    obj = copy.deepcopy(obj)
    print(f"Minimum length for filtering out disconnected HPC nodes : {min_hpc_len}bp")
    
    node = obj.node.copy()
    edge = obj.edge.copy()
    node.reset_index(drop=True, inplace=True)
    edge.reset_index(drop=True, inplace=True)

    nodelst = set(node['node'].tolist())
    edgelst = set(edge['node1'].tolist() + edge['node2'].tolist() )

    disconnected_nodes = list(nodelst - edgelst)
    print(f"Number of disconnected nodes found : {len(disconnected_nodes)}")

    disconnected_nodes_db = node.loc[node['node'].isin(disconnected_nodes),:]
    disconnected_nodes_db_rm = disconnected_nodes_db.loc[disconnected_nodes_db['len'] < min_hpc_len, :]
    disconnected_nodes_db_rm = disconnected_nodes_db_rm['node'].tolist()
    print(f"Number of disconnected nodes that will be filtered : {len(disconnected_nodes_db_rm)}")

    disconnected_nodes_db_keep = node.loc[node['node'].isin(disconnected_nodes),:]
    disconnected_nodes_db_keep = disconnected_nodes_db_keep.loc[disconnected_nodes_db_keep['len'] >= min_hpc_len, :]
    disconnected_nodes_db_keep = disconnected_nodes_db_keep['node'].tolist()
    print(f"Disconnected paths that exceed a certain length will be preserved : {len(disconnected_nodes_db_keep)}")

    path = obj.paths.copy()
    if "rm" not in path.columns:
        path['rm'] = ""

    path['tmp'] = list(path['path'].str.rstrip("+-"))
    path.loc[path['tmp'].isin(disconnected_nodes_db_rm), 'rm'] = "disconnected_node"
    path.loc[path['tmp'].isin(disconnected_nodes_db_keep), 'rm'] = "keep_long_disconnected_node"
    del path['tmp']
    obj.paths = path
    return obj

def keepContig(obj, contig_lst = None, path_list = None):
    """
    This function keeps specific contigs in the paths DataFrame and marks others for removal.

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'paths' DataFrame in obj.paths.
    contig_lst
        A list of contig names to keep. Default is None.
    path_list
        A list of path names to keep. Default is None.
    Returns
    -------
        The updated verkko fillet object with specific contigs kept in the paths DataFrame.
    """
    
    obj = copy.deepcopy(obj)
    path = obj.paths.copy()
    scfmap = obj.scfmap.copy()
    if "rm" not in path.columns:
        path['rm'] = ""
    
    if contig_lst is not None:
        path_name = scfmap[scfmap['contig'].isin(contig_lst)]['pathName'].to_list()
    
    if path_name is not None and path_list is not None:
        path_name = path_name + path_list
    
    path.loc[path['name'].isin(path_name),'rm'] = "keep_contig"

    obj.paths = path
    return obj


def updateConnect(obj):
    """
    This function updates the paths DataFrame by adding a 'rm' column to mark specific paths for removal.

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'paths' DataFrame in obj.paths.
    Returns
    -------
        The updated verkko fillet object with the 'rm' column added to the paths DataFrame.
    """
    obj = copy.deepcopy(obj)
    path = obj.paths.copy()
    gaps = obj.gaps.copy()

    gaps = gaps.loc[gaps['gaps'].isin(['startMarker', 'endMarker'])].reset_index(drop=True)

    for i in tqdm(range(len(gaps))):
        if gaps['cat'][i] == "connectContig":
            note = gaps.loc[i, "notes"].split(" ")
            source = note[1]
            target = note[3]
            ori = note[5]
            flip = note[8]
            
            path.loc[path['name'] == source, 'rm'] = 'rm-connect-source' + "_" + target + "_" + ori + "_" + flip

        elif gaps['cat'][i] == "addContig":
            print("")
    

    obj.paths = path
    return obj

def writeFixedGraph(path, graph, out = "assembly.fixed.paths.gfa", dupNode=[]):
    """
    This function writes the fixed graph to a GFA file.

    Parameters
    ----------
    path
        The path DataFrame containing the paths to be written to the GFA file.
    graph
        The graph file to be used for writing the fixed graph.
    dupNode
        A list of duplicate nodes to be excluded from the graph. Default is an empty list.
    Returns
        None
    """
    print(f"Writing fixed graph to {out}...")
    path = path.copy()
    onlyused = path.loc[path['rm'].str.startswith("keep"), 'path'].to_list()
    onlyused = [item for s in onlyused for item in s.split(',')]
    onlyused = [item.rstrip(r'+-') for item in onlyused]
    # onlyused

    gapNode = list(set(path.loc[path['rm'] == "keep_Nodes_in_unresolved_gaps", 'path'].unique()))
    gapNode = [item.rstrip(r'+-') for item in gapNode]

    node = pd.read_csv(graph, sep='\t', header=None)
    node = node[node[0] == "S"]
    node.columns = ['type', 'node','seq', 'ln', 'rc','il']
    node = node.loc[node['node'].isin(onlyused),:]
    node = node.loc[~node['node'].isin(dupNode),:]
    node.reset_index(drop=True, inplace=True)

    edge = pd.read_csv(graph, sep='\t', header=None)
    edge.reset_index(drop=True, inplace=True)
    edge = edge[edge[0] == "L"]
    edge.columns = ['type', 'node1','node1_strand', 'node2', 'node2_strand','overlap']
    edge = edge.loc[edge['node1'].isin(onlyused),:]
    edge = edge.loc[edge['node2'].isin(onlyused),:]
    edge.reset_index(drop=True, inplace=True)

    contigdf  = path.loc[path['rm'] == "keep_contig",:]
    contigdf.reset_index(drop=True, inplace=True)

    contig = contigdf['path'].to_list()
    contig_name = contigdf['name'].to_list()

    edge['keep'] = False
    dupNodeDf = pd.DataFrame()
    dupEdgeDf = pd.DataFrame()
    # for path in contig:
    for i in tqdm(range(len(contig))):
        path_contig = contig[i].split(",")
        # path_contig = [item.rstrip(r'+-') for item in path_contig]
        contig_name_list = contig_name[i]

        for j in range(len(path_contig)-1):
            # extract last str ( + or -)
            source = path_contig[j]
            source_add = path_contig[j]
            source_strand = "-"
            if source.endswith("-"):
                source_strand = "-"
            source = source.rstrip(r"+-")
            source_add = source_add.rstrip(r"+-")

            target = path_contig[j+1]
            target_add = path_contig[j+1]
            target_strand = "_"
            if target.endswith("-"):
                target_strand = "-"
            target = target.rstrip(r"+-")
            target_add = target_add.rstrip(r"+-")
            
            # make unique node
            if source in dupNode:
                source_add = source + "_" + contig_name_list
                
                subDf = node.loc[node['node'] == source]
                subDf['node'] = subDf['node'].str.replace(source, source_add)
                dupNodeDf = pd.concat([dupNodeDf, subDf], ignore_index=True)

            if target in dupNode:
                target_add = target + "_" + contig_name_list
                
                subDf = node.loc[node['node'] == target]
                subDf['node'] = subDf['node'].str.replace(target, target_add)
                dupNodeDf = pd.concat([dupNodeDf, subDf], ignore_index=True)

            if source in dupNode or target in dupNode:
                #print(source, target)
                subedge = edge.loc[(edge['node1'] == source) & (edge['node2'] == target)]
                #print(subedge)
                subedge['node1'] = subedge['node1'].str.replace(source, source_add)
                subedge['node2'] = subedge['node2'].str.replace(target, target_add)
                #print(subedge)
                edge = pd.concat([edge, subedge], ignore_index=True)

            if source.startswith("[") or target.startswith("["):
                edge.loc[edge['node1'].isin([source_add, target_add]) | edge['node2'].isin([source_add, target_add]), 'keep'] = True
            else:
                edge.loc[edge['node1'].isin([source_add, target_add]) & edge['node2'].isin([source_add, target_add]), 'keep'] = True
                
                # print(source, target)

    edge.loc[edge['node1'].isin(gapNode) | edge['node2'].isin(gapNode), 'keep'] = True

    edge = edge[edge['keep'] == True]
    edge.drop_duplicates(inplace=True)
    del edge['keep']
    edge = edge.loc[~edge['node1'].isin(dupNode),:]
    edge = edge.loc[~edge['node2'].isin(dupNode),:]
    edge.reset_index(drop=True, inplace=True)

    node = pd.concat([node, dupNodeDf], ignore_index=True)
    node = node.loc[~node['node'].isin(dupNode),:]
    node.drop_duplicates(inplace=True)
    node.reset_index(drop=True, inplace=True)

    node.to_csv(out, header = False, index = False, sep = '\t')
    edge.to_csv(out, header = False, index = False, sep = '\t', mode='a')