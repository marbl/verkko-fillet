import subprocess
import pandas as pd
import plotly.graph_objects as go
import os
import numpy as np
import itertools
import copy
import pandas as pd
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from sklearn.linear_model import LinearRegression
import math

def impute_depth(obj, loopNode, realNum):
    
    if len(loopNode) != len(realNum):
        print(f"loopNum : {len(loopNode)}")
        print(f"realNum : {len(realNum)})")
        raise ValueError("The length of loopNode and realNum should be the same")
    depthDb= obj.node.copy()
    depthDb.index= depthDb['node']
    depthdb_sub = depthDb.loc[loopNode]
    depthdb_sub

    # Define the lists
    list1 = list(depthdb_sub['norm_len'])

    list2 = realNum

    # Remove NaNs from list2 and corresponding elements from list1 for training
    train_indices = [i for i, x in enumerate(list2) if not math.isnan(x)]
    train_list1 = [list1[i] for i in train_indices]
    train_list2 = [list2[i] for i in train_indices]

    # Reshape data for linear regression
    train_list1 = np.array(train_list1).reshape(-1, 1)
    train_list2 = np.array(train_list2)

    # Fit linear regression model
    model = LinearRegression()
    model.fit(train_list1, train_list2)

    # Impute the NaN values
    imputed_values = []
    for i, x in enumerate(list2):
        if math.isnan(x):
            # Predict using the model
            imputed_value = model.predict(np.array([[list1[i]]]))[0]
            imputed_values.append(imputed_value)
        else:
            imputed_values.append(x)

    # print(imputed_values)

    depthdb_sub['imputed'] = imputed_values
    depthdb_sub['real'] = realNum
    return depthdb_sub

def calNodeDepth(obj, width = 8, height = 5, save = True, figName = None, dpi = 300, force = False):
    """
    Plot the depth of nodes in the graph.

    Parameters
    -----------
    obj : object
        The VerkkoFillet object.
    width : int
        The width of the plot. Default is 8.
    height : int
        The height of the plot. Default is 5.
    save : bool
        Whether to save the plot. Default is True.
    figName : str
        The name of the figure file. Default is None.
    dpi : int
        The resolution of the figure. Default is 300.
    """

    
    obj = copy.deepcopy(obj)

    if obj.paths_freq is None:
        raise ValueError("No paths_freq found. Please run `vf.tl.graphIdx` and `vf.tl.graphAlign` first.")
    if obj.node is None:
        raise ValueError("No node information found. Please run `vf.tl.graphIdx` and `vf.tl.graphAlign` first")
    columnList= ["count","norm_len","cov","cov_hap"]
    if all(elem in obj.node.columns for elem in columnList):
        if force:
            for col in columnList:
                del obj.node[col]
            else:
                print("The depth of nodes has already been calculated. Skip calculation.")
        
        return obj
    

    # Count the frequency of each node in the paths
    print("Counting the frequency of each node in the paths...")
    paths_freq = obj.paths_freq.copy()
    
    # Ensure `path_modi` is a string before splitting
    paths_freq['nested_list'] = paths_freq['path_modi'].astype(str).str.split("@")

    # Fix the apply function
    paths_freq['expanded_list'] = paths_freq.apply(
        lambda row: [item for item in row['nested_list'] for _ in range(row['nsupport'])], axis=1
    )

    # Flatten the list while filtering out None values
    flatlist = list(filter(None, itertools.chain(*paths_freq['expanded_list'])))

    # Count occurrences
    count_dict = Counter(flatlist)
    count_df = pd.DataFrame(count_dict.items(), columns=['node', 'count'])

    # Sort and get median
    count_df = count_df.sort_values(by='count', ascending=False)

    nodedb = obj.node.copy()

    color_map = {
        "#AAAAAA" : "ambiguous",
        "#FF8888" : "hap1",
        "#8888FF" : "hap2",
        "#FFFF00" : "ambiguous"
    }
    nodedb['hap'] = nodedb['color'].map(color_map)
    # hapdb

    merged = nodedb[['node', 'len', 'hap']]
    merged = pd.merge(merged, count_df, on='node')
    merged['norm_len'] = merged['count'] / merged['len']
    
    median = merged['norm_len'].median()
    merged['cov'] = merged['norm_len']/median
    print(f"Median normalized by length: {median}")

    hapSpecificMed = merged.loc[merged['hap'].isin(["hap1", "hap2"]),"norm_len"].median()
    merged['cov_hap'] = merged['norm_len'] / hapSpecificMed
    print(f"Median of haplotype specific nodes normalized by length : {hapSpecificMed}")
    
    count_df = merged.copy()
    y = 'count'
    count_df = count_df.sort_values(by=y, ascending=False)

    # plot line plots
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1])  # Line plot wider than box plot

    # Line plot
    print(f"Plotting the depth of nodes in the graph...")
    
    ax0 = fig.add_subplot(gs[0])
    sns.lineplot(data=count_df, x='node', y=y, marker="o", ax=ax0, markersize=3, alpha=0.5)
    ax0.set_title("Line Plot")
    ax0.set_yscale("log")
    ax0.set(xticklabels=[], xlabel=None, xticks=[])
    ax0.axhline(count_df[y].median(), color='r', linestyle='--', label='Median', linewidth=.5, alpha=0.5)

    # Box plot (on the right)
    ax1 = fig.add_subplot(gs[1])
    sns.boxplot(y=count_df[y], ax=ax1)
    ax1.set_title("Box Plot")
    ax1.set_yscale("log")

    plt.title("Depth of nodes in the graph")
    plt.tight_layout()
    
    if figName is None:
        figName = f"figs/nodeDepth.png"

    if save:
        if not os.path.exists("figs"):
            os.makedirs("figs")

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name")

        elif not os.path.exists(figName):
            plt.savefig(figName, dpi=dpi)
            print(f"File {figName} saved")

    plt.show()
    #print(count_df)
    #print(nodedb)
    count_df.reset_index(drop=True, inplace=True)
    nodedb = pd.merge(nodedb, count_df.loc[:,['node'] + columnList], on='node', how='outer')
    obj.node = nodedb

    return obj

def estLoops(obj, nodeList, gaf="graphAlignment/verkko.graphAlign_allONT.gaf"):
    """
    Estimate the number of loops between two nodes in the graph.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    nodeList
        A list of two nodes to be compared.
    gaf
        Path to the GAF file containing the graph alignment information. Default is "graphAlignment/verkko.graphAlign_allONT.gaf".   
    """
    if obj is None:
        raise ValueError("The 'obj' parameter is required.")
    if nodeList is None or len(nodeList) < 2:
        raise ValueError("The 'nodeList' parameter must contain two nodes.")

    
    node1, node2 = nodeList
    gaf = os.path.abspath(gaf)

    if not os.path.exists(gaf):
        raise FileNotFoundError(f"The GAF file does not exist: {gaf}")
    
    # Step 1: Run the first grep command to search for node2 in the file
    grep1 = subprocess.run(['grep', '-w', node2, gaf], capture_output=True, text=True)
    
    if grep1.returncode != 0:  # Check if grep1 failed
        print(f"No reads detected for node2: {node2}")
        return
    
    # Step 2: Run the second grep command to search for node1 in the output of the previous command
    grep2 = subprocess.run(['grep', '-w', node1], input=grep1.stdout, capture_output=True, text=True)
    
    if grep2.returncode != 0:  # Check if grep2 failed
        print(f"No reads detected for node1: {node1}")
        return
    
    # Proceed with the rest of the pipeline if both grep commands succeed
    try:
        sed1 = subprocess.run(['sed', 's/id:f://g'], input=grep2.stdout, capture_output=True, text=True)
        awk1 = subprocess.run(['awk', '{if ($NF > 0.99 && $3 < 100 && $4 + 100 > $2) print $6}'], input=sed1.stdout, capture_output=True, text=True)
        tr1 = subprocess.run(['tr', '>', ','], input=awk1.stdout, capture_output=True, text=True)
        tr2 = subprocess.run(['tr', '<', ','], input=tr1.stdout, capture_output=True, text=True)
        sed2 = subprocess.run(['sed', f's/{node1}//g'], input=tr2.stdout, capture_output=True, text=True)
        sed3 = subprocess.run(['sed', f's/{node2}//g'], input=sed2.stdout, capture_output=True, text=True)
        sed4 = subprocess.run(['sed', 's/,,//g'], input=sed3.stdout, capture_output=True, text=True)
        awk2 = subprocess.run(['awk', '-F', ',', '{print NF-1}'], input=sed4.stdout, capture_output=True, text=True)
        sort1 = subprocess.run(['sort'], input=awk2.stdout, capture_output=True, text=True)
        uniq = subprocess.run(['uniq', '-c'], input=sort1.stdout, capture_output=True, text=True)
        awk3 = subprocess.run(['awk', '{print $2, $1}'], input=uniq.stdout, capture_output=True, text=True)
        final_output = subprocess.run(['sed', '-e', f"1i {node1}_to_{node2}"], input=awk3.stdout, capture_output=True, text=True)
    except Exception as e:
        print(f"Error during processing: {e}")
        return
    
    # Process and plot results
    lines = final_output.stdout.strip().split('\n')
    header = [lines[0]]  # First line contains the header
    rows = [line.split() for line in lines[1:]]  # Remaining lines are data
    
    df = pd.DataFrame(rows, columns=['nLoop', 'freq'])
    df['nLoop'] = pd.to_numeric(df['nLoop'])
    df['freq'] = pd.to_numeric(df['freq'])
    df = df.sort_values(by="nLoop")
    
    all_nLoops = pd.DataFrame({'nLoop': range(df['nLoop'].min(), df['nLoop'].max() + 1)})
    df = pd.merge(all_nLoops, df, on='nLoop', how='left').fillna({'freq': 0})
    
    # Display DataFrame
    # print(df)
    
    # Plotting
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df['nLoop'], 
        y=df['freq'],
        mode='lines+markers',
        marker=dict(symbol='circle', size=10, color='blue'),
        name='Freq'
    ))
    fig.update_layout(
        title=str(header),
        xaxis_title='Number of loops',
        yaxis_title='Frequency',
        showlegend=False,
        width=500,
        height=400,
        plot_bgcolor='white',
        xaxis=dict(showgrid=True, gridcolor='grey'),
        yaxis=dict(showgrid=True, gridcolor='grey')
    )
    fig.show()

def checkNodeInfo(obj, node_list):
    """
    Check if the nodes in the node_list are in the obj.node['node'] and print the information if they are.
    Parameters
    ----------
    obj : VerkkoFillet object
        The object containing the node information.
    node_list : list
        A list of nodes to be checked.
    Returns
    -------
    node_info : DataFrame
        A DataFrame containing the information of the nodes in the node_list that are in the obj.node['node'].
    """
    # check if node_name is in obj.node['node'] and print the information if it is
    node_in = []
    node_out = []
    for node in node_list:
            if node in obj.node['node'].values:
                node_in.append(node)
            else:
                node_out.append(node)
    
    print(f"Nodes in obj.node: {node_in}")
    if len(node_out) > 0:
        print(f"Nodes not in obj.node: {node_out}")
    
    node_info = obj.node[obj.node['node'].isin(node_in)]
    
    return node_info