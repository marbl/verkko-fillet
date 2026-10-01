import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from .._run_shell import run_shell
import matplotlib.cm as cm
import copy
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def showMashmapOri(obj, mashmap_out = "chromosome_assignment/assembly.mashmap.out", by = "chr_hap", 
                   width = 5, height = 7, save = True, figName = None, force = False):
    """
    Generates a bar plot showing the covered regions of the assembly for each reference.

    Parameters
    -----------
    obj
        An object that contains a .stats attribute, which should be a pandas DataFrame.
    mashmap_out
        The Mashmap output file, aligning the assembly to the reference. Default is assembly.mashmap.out.filtered.out.
    by
        Specifies the grouping method for the plot. Default is "chr_hap". Available options are ['contig', 'all', 'chr_hap'].
    """
    obj = copy.deepcopy(obj)
    working_dir = os.path.abspath(obj.verkko_fillet_dir)  # Ensure absolute path for the working directory

    print(f"Reading Mashmap output from {working_dir + '/' + mashmap_out}")
    mashmap = pd.read_csv(working_dir + "/" + mashmap_out , header = None, sep ='\t')
    mashmap.columns = ['qname','qlen','qstart','qend','strand','tname','tlen','tstart','tend','nmatch','blocklen','mapQ','id','kc']
    mashmap['block_q'] =  mashmap['qend'] - mashmap['qstart']
    grouped = mashmap.groupby(['qname', 'tname', 'strand'])
    data = grouped.agg(
        qlen=('qlen', 'first'),  # Take the first value of qlen as representative
        qcover=('block_q','sum')  # Calculate coverage
    ).reset_index()
    
    # Calculate percentage coverage
    data['qcover_perc'] = data['qcover'] / data['qlen'] * 100
    # Copy the stats DataFrame from obj
    stats = obj.stats.copy()
    # Filter rows based on 'qname' matching contig names in obj.stats['contig']
    contig_list = list(stats['contig'])  # Assuming this is a list of contig names
    data = data.loc[data['qname'].str.contains('|'.join(contig_list)), :]

    # Create a new column 'name' by concatenating 'ref_chr' and 'hap'
    if by == 'chr_hap':
        stats['by'] = stats['ref_chr'].astype(str) + "_" + stats['hap'].astype(str)
    if by == 'contig':
        stats['by'] = stats['contig'].astype(str)
    if by == 'all':
        stats['by'] = stats['contig'].astype(str) + '_' + stats['ref_chr'].astype(str) + "_" + stats['hap'].astype(str)
        
    # Display the first two rows of the DataFrame
    data = pd.merge(data,stats,how= 'left',left_on = 'qname', right_on = 'contig')
    # Drop rows with no matching stats entry, since 'by' would be NaN (float) and break barh
    if data['by'].isna().any():
        missing = data.loc[data['by'].isna(), 'qname'].unique()
        print(f"Warning: no matching stats entry for qname(s): {list(missing)}, dropping from plot")
        data = data.dropna(subset=['by'])
    data['by'] = data['by'].astype(str)
    data.loc[data['strand'] == "-", 'qcover_perc'] *= -1
    # Separate positive and negative values for clarity
    data['positive_qcover'] = data['qcover_perc'].where(data['strand'] == '+', 0)
    data['negative_qcover'] = data['qcover_perc'].where(data['strand'] == '-', 0)
    
    # Sort the data by qname for better visualization
    data = data.sort_values(by = "qcover_perc", ascending = False)
    
    # Prepare the plot
    fig, ax = plt.subplots(figsize=(width, height))
    print(data.head())
    # Add horizontal bars for positive and negative values
    ax.barh(data['by'], data['positive_qcover'], color='purple', label='Positive Strand')
    ax.barh(data['by'], data['negative_qcover'], color='skyblue', label='Negative Strand')
    
    # Add labels, legend, and gridlines
    ax.set_xlabel('qcover_perc (%)')
    ax.set_ylabel('Contig')
    ax.set_title('Horizontal Stacked Bar Plot with Positive and Negative Strands')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for reference
    ax.legend()
    
    # Adjust layout and show the plot
    plt.tight_layout()
    if figName is None:
        figName = f"figs/showMashmapOri.heatmap.pdf"

    if save:
        if not os.path.exists("figs"):
            print("Creating figs directory")
            os.makedirs("figs")

        if force and os.path.exists(figName):
            # remove original plot
            print(f"File {figName} already exists and force is set to True, removing the original file")
            os.remove(figName)

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name")

        elif not os.path.exists(figName):
            plt.savefig(figName, bbox_inches='tight', pad_inches=0.1)
            print(f"File {figName} saved")
    plt.show()

def nodeMashmapBlockSize(
        mashmap_out,
        node,
        showNum = 10,
        width = 8, height = 5,
        save = True, figName = None, force = False):
    """
    Plot the top 10 nodes with the largest block size aligned to a specific node.

    Parameters
    -----------
    mashmap_out : str
        The Mashmap output file.
    node : str
        The node to which the other nodes are aligned.
    showNum : int
        The number of nodes to display. Default is 10.
    width : int
        The width of the plot. Default is 8.
    height : int
        The height of the plot. Default is 5.
    save : bool
        Whether to save the plot. Default is True.
    figName : str
        The name of the figure file. Default is None.
    """
    tab = pd.read_csv(mashmap_out, sep="\t", header=None, usecols=[0,4,5,10])
    tab.columns = ['query', 'strand', 'ref', 'blocksize']
    tab['ref'] = tab['strand'] + tab['ref']
    tab = tab.loc[tab['query'] == node]
    tab_group = tab.groupby('ref')['blocksize'].sum().sort_values(ascending=False).reset_index()
    tab_group.columns = ['node', 'blocksize']
    tab_group = tab_group.head(showNum)

    plt.figure(figsize=(width, height))  # Adjust figure size

    ax = sns.barplot(x="node", y="blocksize", data=tab_group, color="grey")

    # Labels and title
    plt.xlabel("nodes")
    plt.ylabel("block size")
    plt.title(f"Top 10 nodes with the largest block size aligned to {node}")
    plt.xticks(rotation=75)
    # Show the plot
    if figName is None:
        figName = f"figs/mashmap.{node}_top{showNum}_blocksize.pdf"

    if save:
        if not os.path.exists("figs"):
            os.makedirs("figs")

        if force and os.path.exists(figName):
            # remove original plot
            print(f"File {figName} already exists and force is set to True, removing the original file")
            os.remove(figName)

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name")

        elif not os.path.exists(figName):
            plt.savefig(figName, bbox_inches='tight', pad_inches=0.1)
            print(f"File {figName} saved")
    plt.show()