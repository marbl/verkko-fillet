import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import copy


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def plotHist_readOnNode(nodeinfo, mergedb_all, loc_on_node, node, width = 5, height = 7, save = True, figName = None, **kwargs):
    """
    Generates a histogram showing the distribution of reads on a node.

    Parameters
    -----------
    nodeinfo
        A DataFrame containing information about the nodes.
    mergedb_all
        A DataFrame containing the merged data.
    loc_on_node
        A DataFrame containing the location of the node.
    node
        The node to plot.
    width
        Width of the plot. Default is 5.
    height
        Height of the plot. Default is 7.
    save
        If True, the plot is saved as a PNG file. Default is True.
    figName
        Name of the saved plot. Default is None. If None, the plot is saved as "figs/intra_telo.heatmap.png".
    kwargs
        Additional arguments to pass to seaborn.histplot.
    """
    
    mergedb = mergedb_all.loc[mergedb_all['node'] == node]
    hap = nodeinfo.loc[nodeinfo['node'] == node, 'hap_node'].values[0]
    hap_ratio = nodeinfo.loc[nodeinfo['node'] == node, 'mat:pat'].values[0]

    sns.histplot(data = mergedb, bins = 100, kde =True, **kwargs).set_title(f"""{node} ({hap} {hap_ratio})""")

    start_coor = loc_on_node.loc[loc_on_node['node'] == node, 'start'].values[0] 
    end_coor = loc_on_node.loc[loc_on_node['node'] == node, 'end'].values[0]

    plt.figure(figsize=(width, height))
    plt.axvline(x=start_coor , color='red', linestyle='--')
    plt.axvline(x=end_coor , color='blue', linestyle='--')
    if figName is None:
        figName = f"figs/intra_telo.heatmap.png"

    if save:
        if not os.path.exists("figs"):
            print("Creating figs directory")
            os.makedirs("figs")

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name")

        elif not os.path.exists(figName):
            plt.savefig(figName)
            print(f"File {figName} saved")
    plt.show()