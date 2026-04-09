import re
import networkx as nx
from collections import Counter
import pandas as pd
import copy
import sys
import shlex
import subprocess
import os
import shutil
from tqdm import tqdm
from .._default_func import flatten_and_remove_none

def make_unique(column):
    counts = {}
    result = []
    for name in column:
        if name in counts:
            counts[name] += 1
            result.append(f"{name}_{counts[name]}")
        else:
            counts[name] = 0
            result.append(name)
    return result
    
def create_pairs(input_str):
    # Split the string into a list
    split_values = input_str.split(',')
    # Generate pairs of consecutive elements
    return [(split_values[i], split_values[i+1]) for i in range(len(split_values)-1)]

def remove_elements_starting_with_bracket(lst):
    return [item for item in lst if not item.startswith('[')]

def remove_ignore_nodes(lst, ignore_lst):
    return [item for item in lst if item not in ignore_lst]

def find_multi_used_node(obj):
    """\
    Find nodes that are used in more than one path.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    
    Returns
    -------
    duplication reads : list 
        list of duplicated nodes
    """
    if obj.paths is None:
        raise ValueError("Paths are not available in the VerkkoFillet object.")
    if obj.scfmap is None:
        raise ValueError("SCF map is not available in the VerkkoFillet object.")
    if obj.stats is None:
        raise ValueError("Stats are not available in the VerkkoFillet object.")
    path = obj.paths.copy()
    scfmap = obj.scfmap.copy()
    stats = obj.stats.copy()
    
    path['nodeSet'] = path['path'].apply(
        lambda x: [y for y in re.split(r'[><]', x) if y] if re.search(r'[><]', x) else [word.rstrip('-+') for word in x.replace(',', ' ').split() if word.rstrip('-+')]
    )
    path = pd.merge(scfmap, path, how = 'outer',left_on = "pathName",right_on="name")
    path['nodeSet'] = path['nodeSet'].apply(remove_elements_starting_with_bracket)
    
    assignedContig = list(stats['contig']) 
    unassignedPath = path.loc[~path['contig'].isin(assignedContig)].reset_index()
    assignedPath = path.loc[path['contig'].isin(assignedContig)].reset_index()
    
    assignedContig = list(stats['contig'])
    assignedPath = pd.merge(stats, assignedPath, on = 'contig', how = 'outer')
    
    path_grouped = assignedPath.groupby('ref_chr')['nodeSet'].agg(
        lambda x: set([item for sublist in x if isinstance(sublist, (list, set)) for item in sublist])
    ).reset_index()
    
    list_of_lists = path_grouped['nodeSet'].tolist()
    
    # Example list of lists
    
    # Flatten the list of lists and count the occurrences of each element
    flat_list = [item for sublist in list_of_lists for item in sublist]
    element_counts = Counter(flat_list)
    
    # Find the elements that appear in more than one list
    duplicates = [item for item, count in element_counts.items() if count > 1]
    
    return duplicates, path_grouped

def naming_contigs(obj, node_database, duplicate_nodes ,
                   dam = "mat", sire = "pat", fai = "assembly.fasta.fai"):
    """\
    Rename the contigs based on the provided chromosome map file.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    node_database
        The DataFrame containing the mapping of nodes to chromosomes.
    duplicate_nodes
        List of duplicated nodes.
    gfa
        The path to the GFA file. Default is "assembly.homopolymer-compressed.noseq.gfa".
    dam
        The name of the dam. Default is "mat".
    sire
        The name of the sire. Default is "pat".
    fai
        The path to the FASTA index file. Default is "assembly.fasta.fai".
    
    Returns
    -------
        The DataFrame containing the nodes and their corresponding assigned contig names.

    """
    obj = copy.deepcopy(obj)
    if obj.stats is None:
        raise ValueError("Stats are not available in the VerkkoFillet object.")
    if obj.paths is None:
        raise ValueError("Paths are not available in the VerkkoFillet object.")
    if obj.scfmap is None:
        raise ValueError("SCF map is not available in the VerkkoFillet object.")
    if obj.node is None:
        raise ValueError("Nodes are not available in the VerkkoFillet object.")
    if obj.edge is None:
        raise ValueError("Edges are not available in the VerkkoFillet object.")
    
    stats = obj.stats.copy()
    path = obj.paths.copy()
    scfmap = obj.scfmap.copy()
    gfa_link = obj.edge.copy()
    gfa_link = gfa_link[['node1', 'node2']]
    gfa_link.columns = ['start', 'end']

    print(f"All required files are available ... ")
    
    if not os.path.exists(fai):
        raise FileNotFoundError(f"The file {fai} does not exist.")
    else:
        print(f"Reading fai file ... ")
        fai_chr = pd.read_csv(fai, sep='\t', header=None, usecols=[0])[0].tolist()

    scfmap = scfmap.loc[scfmap['contig'].isin(fai_chr)]
    
    gfa_link = gfa_link[~gfa_link['start'].isin(duplicate_nodes)]
    gfa_link = gfa_link[~gfa_link['end'].isin(duplicate_nodes)]
    
    path['nodeSet'] = path['path'].apply(
        lambda x: set(y for y in re.split(r'[><]', x) if y) if re.search(r'[><]', x) else set(word.rstrip('-+') for word in x.replace(',', ' ').split() if word.rstrip('-+'))
    )
    path['path'] = path['nodeSet'].apply(lambda x: ','.join(x))
    path = pd.merge(scfmap, path, how = 'left',left_on = "pathName",right_on="name")
    path['nodeSet'] = path['nodeSet'].apply(remove_elements_starting_with_bracket)
    path['nodeSet'] = path['nodeSet'].apply(lambda lst: remove_ignore_nodes(lst, duplicate_nodes))
    
    assignedContig = list(stats['contig'])
    assignedContig = [re.sub(r':.*', '', string) for string in assignedContig]
    unassignedPath = path.loc[~path['contig'].isin(assignedContig)].reset_index()
    assignedPath = path.loc[path['contig'].isin(assignedContig)].reset_index()
    
    # Exploding the list into separate rows
    df_exploded = node_database.explode('nodeSet').reset_index(drop=True)
    
    # Renaming columns to match the desired output
    df_exploded.columns = ['chr', 'node']

    result_df = gfa_link.copy()

    # Create a graph
    G = nx.Graph()
    G.add_edges_from(result_df.values)
    
    # Find connected components
    connected_components = [set(component) for component in nx.connected_components(G)]
    num_cluster = len(connected_components)
    print(f"Number of connected components: {num_cluster}")
    
    # assign chromosome
    dict1 = {}
    
    for i in tqdm(range(0,len(connected_components)), desc="Assigning Chromosome to clusters", ncols=80, colour="white"):
        dict2 = {}
        chr_assign = df_exploded.loc[df_exploded['node'].isin(connected_components[i]), "chr"].unique()
        if len(chr_assign) == 1:
            dict2 = {chr_assign[0] : connected_components[i]}
            dict1.update(dict2)
            #print("component_" + str(i) + " : " + chr_assign.astype(str))
        if len(chr_assign) >1: 
            chrName = "_".join(chr_assign)
            dict2 = {chrName : connected_components[i]}
            dict1.update(dict2)
            #print("component_" + str(i) + " : " + chr_assign)
        if len(chr_assign) < 1:
            print("NodeCluster_" + str(i) + " : empty")
    
    # Assuming dict1 has sets or lists as values and we want to check if 'nodeSet' is a subset of any of those sets
    for i in range(0, unassignedPath.shape[0]):
        # Access the nodeSet for the current row and convert it to a set
        node_set = set(unassignedPath.loc[i, "nodeSet"])
        
        # Find the key that corresponds to a matching value in dict1
        some_key = None  # Initialize the key
        
        for key, value in dict1.items():
            if node_set.issubset(set(value)):  # Check if node_set is a subset of value
                some_key = key  # If a match is found, assign the key
                break  # No need to continue checking after the first match
        
        # If a match was found, assign the key to the 'assignChr' column
        # if some_key is not None:
        unassignedPath.loc[i, "assignChr"] = some_key
        
        # Print if it's a subset and which key was assigned
        # print(f"Row {i}: Is subset? {some_key is not None}, Assigned Key: {some_key}")
    
    print(f"Starting Naming contigs ...")
    # update unassigned
    unassignedPath.loc[unassignedPath['assignChr'].isna(), 'assignChr'] = "chrUn"
    unassignedPath = unassignedPath[~unassignedPath['contig'].isna()].reset_index()
    del unassignedPath['index']
    
    # Split the 'contig' column and extract the first part
    unassignedPath['hap'] = unassignedPath['contig'].str.split(pat="_", expand=True)[0]
    
    # Find rows where 'hap' starts with 'unassigned'
    unassigned_rows = unassignedPath['hap'].str.startswith('unassigned')
    
    # Update the 'hap' column for these rows
    unassignedPath.loc[unassigned_rows, 'hap'] = "hapUn"
    unassignedPath['hap'] = (
        unassignedPath['hap']
        .str.replace('sire', sire, regex=False)
        .str.replace('dam', dam, regex=False)
    )
    unassignedPath['path_id'] = unassignedPath['name'].apply(lambda x: re.sub(r".*_utig", "utig", x))
    unassignedPath['new_contig_name'] = (
        unassignedPath['assignChr'] + "_" + unassignedPath['hap'] + "_random_" + unassignedPath['path_id']
    )

    assignedPath = pd.merge(assignedPath,stats, on='contig', how = 'left')
    assignedPath['new_contig_name'] = assignedPath['ref_chr'].astype(str) + "_" + assignedPath['hap'].astype(str)
    assignedPath= assignedPath[['contig','new_contig_name']]

    unassignedPath= unassignedPath[['contig','new_contig_name']]

    final_contigNaming = pd.concat([assignedPath, unassignedPath])
    final_contigNaming['new_contig_name'] = make_unique(final_contigNaming['new_contig_name'])
    # final_contigNaming.to_csv(out_mapFile, sep ='\t', header = None, index=False)
    # Display the updated DataFrame
    print(f"Done!")
    return final_contigNaming

def cut_graph_using_ancestors(graph, source, target):
    """
    Cuts the graph based on ancestors of the target and descendants of the source.
    Keeps only the nodes and edges that are part of the path from source to target.

    Parameters:
        graph (nx.DiGraph): The original directed graph.
        source (str): The source node.
        target (str): The target node.

    Returns:
        nx.DiGraph: A new graph containing only the nodes and edges involved in paths from source to target.
    """
    if source not in graph or target not in graph:
        print(f"Source or target node does not exist in the graph.")
        print(f"Add edge from {source} to {target} to the graph.")
        graph.add_edge(source, target)

    
    # Find the ancestors of the target (all nodes that can reach the target)
    ancestors_of_target = nx.ancestors(graph, target)
    ancestors_of_target.add(target)  # Include the target node itself

    # Find the descendants of the source (all nodes that can be reached from the source)
    descendants_of_source = nx.descendants(graph, source)
    descendants_of_source.add(source)  # Include the source node itself

    # Find the intersection of ancestors and descendants
    relevant_nodes = ancestors_of_target.intersection(descendants_of_source)

    # Create a subgraph with only the relevant nodes
    subgraph = graph.subgraph(relevant_nodes).copy()

    return subgraph

def grabNodesInGap(obj, source, target):
    """
    Find nodes in the gap between source and target.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    source
        The source node.
    target
        The target node.

    Returns
    -------
    lst
        List of nodes in the gap.
    """
    # Sample data for edges
    edge = obj.edge.copy()

    # Adding start and end columns as per your structure
    edge['start'] = edge['node1'].astype(str) + edge['node1_strand'].astype(str)
    edge['end'] = edge['node2'].astype(str) + edge['node2_strand'].astype(str)
    edge = edge[['start', 'end']]

    # Create a directed graph
    G = nx.DiGraph()

    # Add edges to the graph
    for _, row in edge.iterrows():
        G.add_edge(row['start'], row['end'])
    
    # Cut the graph using ancestors of target and descendants of source
    cut_graph = cut_graph_using_ancestors(G, source, target)
    if cut_graph is None:
        return 
    else:
        lst = list(cut_graph.nodes)
        return lst

def keepNodesInUnresolvedGaps(obj):
    """\
    Find nodes that are used in more than one path.
    
    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    
    Returns
    -------
    obj
        The updated VerkkoFillet object with new paths dataframe.
    result
        List of nodes kept in unresolved gaps.
    """
    gaps = obj.gaps[['gaps','fixedPath','notes','cat','name']].copy()
    gaps = gaps.reset_index()
    path = obj.paths.copy()
    path = path.reset_index()

    if "rm" not in path.columns:
        path['rm'] = ""
    
    keepGapNode = []

    for i in tqdm(range(len(gaps))):
        if gaps['gaps'][i] == 'startMarker':
            if gaps['cat'][i] == "connectContig":
                note = gaps.loc[i, "notes"].split(" ")
                source = note[1]
                source = path.loc[path['name'] == source, 'path'].values[0].split(",")[-1].replace(" ", "")
                target = note[3]
                target = path.loc[path['name'] == target, 'path'].values[0].split(",")[0].replace(" ", "")
                gapinfo = [source, "[gap]", target]
                # print(f"source: {source}, target: {target}")
            elif gaps['cat'][i] == 'addContig':
                target = gaps['name'][i]
                target_path = path.loc[path['name'] == target, 'path'].values[0]
                fixed_path = gaps['fixedPath'][i]
                gapinfo = [fixed_path, target_path]
        
        elif gaps['gaps'][i] == 'endMarker':
            if gaps['cat'][i] == "connectContig": 
                note = gaps.loc[i, "notes"].split(" ")
                source = note[3]
                source = path.loc[path['name'] == source, 'path'].values[0].split(",")[-1].replace(" ", "")
                target = note[0]
                target = path.loc[path['name'] == target, 'path'].values[0].split(",")[0].replace(" ", "")
                gapinfo = [source, "[gap]", target]
                # print(f"source: {source}, target: {target}")
            elif gaps['cat'][i] == 'addContig':
                target = gaps['name'][i]
                target_path = path.loc[path['name'] == target, 'path'].values[0]
                fixed_path = gaps['fixedPath'][i]
                gapinfo = [target_path, fixed_path]
                print(gapinfo)

        elif gaps['fixedPath'][i] == '':
            gapinfo = gaps['gaps'][i]
        
        else : 
            gapinfo = gaps['fixedPath'][i]
        
        if isinstance(gapinfo, str):
            gapinfo = gapinfo.split(",")

        # print(len(gapinfo))
        # Step 2: Remove all whitespace from each element
        gapinfo = [s.replace(" ", "") for s in gapinfo]

        # Step 3: Find index of first item starting with "["
        index = [i for i, item in enumerate(gapinfo) if item.startswith("[")]
        if index is not None:
            for idx in range(len(index)):
                sliceIndex = index[idx]
                if sliceIndex == 0 or sliceIndex >= len(gapinfo) - 1:
                    continue
                source = gapinfo[sliceIndex-1]
                target = gapinfo[sliceIndex+1]
                # print(f"source: {source}, target: {target}")
                gapNode_to_keep = grabNodesInGap(obj, source, target)
                # gapNode_to_keep= list(gapNode_to_keep)
                keepGapNode.append(gapNode_to_keep)

    result = flatten_and_remove_none(keepGapNode)
    result = [x.rstrip("-+") for x in result]
    print(f"Total number of nodes used in gap filling: {len(result)}")
    print(result)
    includelst = []

    for i in tqdm(range(len(path)), ncols=80, colour="white", desc="Finding paths consisting of unused nodes"):
        path_list = path.loc[i, 'path'].split(",")
        path_list = [x.rstrip("-+") for x in path_list]  # Remove trailing "-" and "+" from each string
        path_name = path.loc[i, 'name']
        # print(path_list)
        if all(item in result for item in path_list):
            includelst.append(path_name)
        
    includelst = list(set(includelst))
    print(f"The total number of paths consisting of nodes should be preserved.: {len(includelst)}")

    path.loc[path['name'].isin(includelst), 'rm'] = "keep_Nodes_in_unresolved_gaps"
    
    if "index" in path.columns:
        del path['index']
    if "level_0" in path.columns:
        del path['level_0']

    obj.paths = path
    return obj, result

from tqdm import tqdm 

def reClusteringGapNodeByPath(obj):
    gapNodeDb =pd.DataFrame()
    paths = obj.paths.copy()
    paths.reset_index(drop=True, inplace=True)

    gaps = obj.gaps.copy()
    gaps.reset_index(drop=True, inplace=True)

    for i in tqdm(range(len(obj.gaps)), desc="Processing gaps", unit="gap"):
        gapinfo = gaps.loc[i, 'gaps']
        gapinfo = [s.replace(" ", "") for s in gapinfo]


        name = gaps.loc[i, 'name']
        hap = name.split("_")[0]
        gapid = gaps.loc[i, 'gapId']
        # print(gapid)

        gap_node = grabNodesInGap(obj, gapinfo[0], gapinfo[-1])
        # print(gap_node)
        if len(gap_node) == 0:
            print(f"Gap node is None for gap ID: {gapid}")
            continue
        gap_node.remove(gapinfo[0])
        gap_node.remove(gapinfo[-1])
        gap_node = [x.rstrip("-+") for x in gap_node]
        

        for j in gap_node:
            node = j + "+"
            # print("Checking node:", node)

            match = paths.loc[paths['path'] == node, 'name']
            
            if not match.empty:
                paths_name = match.values[0]
                # print("Matched path name:", paths_name)

                nodeDf = pd.DataFrame({
                    "gapId": [gapid],
                    "mainContig": [name],
                    "pathName": [paths_name],
                    "node": [node]
                })

                gapNodeDb = pd.concat([gapNodeDb, nodeDf], ignore_index=True)

    scfmap = obj.scfmap.copy()
    scfmap.reset_index(drop=True, inplace=True)
    scfmap = scfmap[['contig','pathName']]
    
    stats = obj.stats.copy()

    stat_scfmap = pd.merge(stats, scfmap, how = 'left', left_on = 'contig', right_on = 'contig')

    # if gapNodeDb is empty, stop
    if gapNodeDb.empty:
        print("The database of gap-nodes is empty.")
        node_feature_dict=None
        stat_scfmap.rename(columns = {'contig' : 'original_contig'}, inplace = True)
        return stat_scfmap,  node_feature_dict
    else:
        gapNodeDb = gapNodeDb.reset_index(drop=True)

        gapNodeDb["mainContig_hap"] = gapNodeDb["mainContig"].str.split("_").str[0]
        gapNodeDb["pathName_hap"] = gapNodeDb["pathName"].str.split("_").str[0]
        gapNodeDb.loc[gapNodeDb['mainContig_hap'] == gapNodeDb['pathName_hap'], 'same'] = True
        gapNodeDb.loc[gapNodeDb['mainContig_hap'] != gapNodeDb['pathName_hap'], 'same'] = False


        gapNodeDb_false = gapNodeDb.loc[gapNodeDb['same'] == False,:]
        gapNodeDb_false.reset_index(drop=True, inplace=True)
        gapNodeDb_false_count = gapNodeDb_false.groupby('node')['mainContig'].nunique().reset_index(name='count').sort_values(by='count', ascending=False)
        unique_nodes = gapNodeDb_false_count.loc[gapNodeDb_false_count['count'] == 1, 'node'].tolist()
        notUnique_nodes = gapNodeDb_false_count.loc[gapNodeDb_false_count['count'] > 1, 'node'].tolist()

        gapNodeDb['fixed'] = gapNodeDb['mainContig']
        gapNodeDb.loc[gapNodeDb['node'].isin(unique_nodes),'fixed'] = gapNodeDb.loc[gapNodeDb['node'].isin(unique_nodes),'mainContig']
        gapNodeDb.loc[gapNodeDb['node'].isin(notUnique_nodes),'fixed'] = "unhap"
        gapNodeDb = gapNodeDb[['mainContig', 'pathName', 'node', 'fixed','same']]
        gapNodeDb.drop_duplicates(inplace=True)
        gapNodeDb.reset_index(drop=True, inplace=True)

        gapNodeDb['fixed'] = gapNodeDb['fixed'] + "_" + gapNodeDb['node']

        stat_scfmap_gapNodedb = pd.merge(stat_scfmap, gapNodeDb, how = 'right', left_on = 'pathName', right_on = 'mainContig')
        stat_scfmap_gapNodedb = stat_scfmap_gapNodedb[['node','pathName_x','pathName_y','ref_chr', 'hap','same']]
        stat_scfmap_gapNodedb.columns = ['node','mainContig_pathName','original_pathName', 'ref_chr', 'hap', 'same']
        stat_scfmap_gapNodedb['ref_chr_pathName'] = stat_scfmap_gapNodedb['ref_chr'].astype(str) + "_" + stat_scfmap_gapNodedb['hap'].astype(str) + "_random_" + stat_scfmap_gapNodedb['node'].astype(str)
        stat_scfmap_gapNodedb = pd.merge(stat_scfmap_gapNodedb, scfmap, how = 'left', left_on = 'original_pathName', right_on = 'pathName')
        stat_scfmap_gapNodedb.rename(columns = {'contig' : 'original_contig'}, inplace = True)
        node_feature_dict = dict(zip(stat_scfmap_gapNodedb["original_contig"], stat_scfmap_gapNodedb["ref_chr_pathName"]))

    
    return stat_scfmap_gapNodedb,node_feature_dict