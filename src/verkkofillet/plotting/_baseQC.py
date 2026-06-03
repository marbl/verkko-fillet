import pandas as pd
import os
import copy
import seaborn as sns
import matplotlib
import matplotlib.patches as mpatches
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
from natsort import natsorted

import logging
logging.getLogger('fontTools').setLevel(logging.WARNING)
logging.getLogger('fontTools.subset').setLevel(logging.WARNING)
logging.getLogger('fontTools.ttLib').setLevel(logging.WARNING)


def n50Plot(obj, width = 8,height = 5 , save = True, figName = None, nprint = None, colName = "ref_chr", dpi = 300, force = False):
    """
    Generates a line plot showing the distribution of contig lengths. The N50 value is indicated by a red dashed line.

    Parameters
    -----------
    obj
        An object that contains a `.stats` attribute, which should be a pandas DataFrame.
    width
        Width of the plot. Default is 8.
    height
        Height of the plot. Default is 5.
    save
        If True, the plot is saved as a PDF file. Default is True.
    figName
        Name of the saved plot. Default is None. If None, the plot is saved as "figs/N50.pdf". You can set the path and format to save the plot. For example, "figs/N50.pdf".
    nprint
        Number of contigs to display. Default is None. If None, all contigs are displayed.
    colName
        Column name to group by in obj.stats. Default is "ref_chr".
    dpi
        Resolution of the saved plot. Default is 300.
    """
    
    obj = copy.deepcopy(obj)

    if obj.stats is None:
        raise AttributeError("Object does not have a 'stats' attribute")
    
    stats = obj.stats.copy()
    if nprint is None:
        nprint = stats.shape[0]

    # check the column name
    if colName not in stats.columns:
        raise ValueError(f"Column name {colName} not found in stats")
    if "hap_verkko" not in stats.columns:
        raise ValueError(f"Column name hap_verkko not found in stats")
    
    stats[colName] = stats[colName].astype(str)
    stats['by'] = stats[colName] + "_" + stats["hap_verkko"]
    stats.sort_values("contig_len", ascending=False, inplace=True)
    stats.reset_index(drop=True, inplace=True)
    total = stats["contig_len"].sum()

    # calculate n50 and l50
    n50 = 0
    l50 = 0

    for i in range(len(stats)):
        if n50 < total/2:
            n50 += stats.loc[i, "contig_len"]
            l50 += 1
        else:
            break

    plt.figure(figsize=(width, height))
    sns.lineplot(data = stats.head(nprint), x = 'by', y = "contig_len")
    plt.xlabel("Contig")
    plt.ylabel("Contig length")
    plt.title(f"Contig length distribution. Top {nprint} contigs")
    plt.axvline(x=l50, color='r', linestyle='--')
    plt.text(l50, stats["contig_len"].max() * 0.1, f"n50={n50:,}", rotation=90, ha = 'left')
    plt.xticks(rotation=75, ha='right')

    if figName is None:
        figName = f"figs/N50.lineplot.pdf"
        
    if save:
        if not os.path.exists("figs"):
            os.makedirs("figs")
        
        if force and os.path.exists(figName):
            # remove original plot
            print(f"File {figName} already exists and force is set to True, removing the original file")
            os.remove(figName)

        if os.path.exists(figName):
            print(f"File {figName} already exists")
            print("Please remove the file or change the name to save the new file")

        elif not os.path.exists(figName):
            plt.savefig(figName, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
            print(f"File {figName} saved")
    print(f"n50={n50:,}")
    print(f"l50={l50}")

    plt.show()


def qvPlot(obj, width = 5, height = 7, save = True, figName = None):
    """
    Generates a bar plot showing QV stats by haplotype.

    Parameters
    -----------
    obj
        An object that contains a `.stats` attribute, which should be a pandas DataFrame.
    """
    # Create figure and axes
    obj = copy.deepcopy(obj)
    qvTab=obj.qv.copy()
    
    fig, ax1 = plt.subplots(figsize=(width, height))
    
    # Create barplot for QV
    barplot = sns.barplot(x='asmName', y='QV', data=qvTab, ax=ax1, color='grey')
    
    # Add labels to the bars
    for index, bar in enumerate(barplot.patches):
        height = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width() / 2,  # X-coordinate (center of bar)
            height-5,                        # Y-coordinate (middle of bar)
            f'{qvTab["QV"][index]:.2f}',       # Label (formatted as an integer)
            ha='center',                       # Horizontal alignment
            va='center',                       # Vertical alignment
            color='white',                     # Text color
            fontsize=10                        # Font size
        )
    
    # Create the second y-axis
    ax2 = ax1.twinx()
    
    # Create a line plot for ErrorRate
    sns.lineplot(x='asmName', y='ErrorRate', data=qvTab, ax=ax2, color='black', label='ErrorRate', marker='o')
    
    # Set axis labels
    ax1.set_ylabel('QV', color='grey')
    ax2.set_ylabel('Error Rate', color='black')
    
    # Set x-axis label
    ax1.set_xlabel('Name of assembly')
    
    # Adjust colors for visibility
    ax1.tick_params(axis='y', colors='grey')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    ax2.tick_params(axis='y', colors='black')
    
    # Remove legend from the second axis (if not needed)
    ax2.legend_.remove() if ax2.legend_ else None
    if figName is None:
        figName = f"figs/qvplot.barplot.pdf"

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
    # Show the plot
    plt.show()

def completePlot(obj, width = 6, height = 4, save = True, figName = None, force = False):
    """
    Generates a bar plot showing contig completeness grouped by reference chromosome and haplotype. The completeness of each chromosome is calculated by comparing it to the reference length. A completeness value greater than 100 indicates that the contig length exceeds the original reference length.

    Parameters
    -----------
    obj
        An object that contains a `.stats` attribute, which should be a pandas DataFrame 
        with the following columns:
        - `ref_chr` (str): Reference chromosome identifier.
        - `hap` (str): Haplotype information.
        - `contig_len` (int): Length of the contigs.
    """
    obj = copy.deepcopy(obj)
    stat_db = obj.stats.copy()
    plt.figure(figsize=(width, height))  # Adjust the figure size as needed
    sns.barplot(stat_db.groupby(['ref_chr','hap'])['completeness'].sum().reset_index(),
                x="ref_chr", y="completeness", hue="hap")
    plt.title("completeness", fontsize=14)
    plt.xticks(rotation=45)
    if figName is None:
        figName = f"figs/completePlot.barplot.pdf"

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


def contigLenPlot(obj, width = 6, height = 3, save = True, figName = None, force = False):
    """
    Generates a bar plot showing length of contig by haplotype.

    Parameters
    -----------
    obj
        An object that contains a `.stats` attribute, which should be a pandas DataFrame.
    """
    obj = copy.deepcopy(obj)
    stat_db = obj.stats.copy()
    plt.figure(figsize=(width, height))  # Adjust the figure size as needed
    sns.barplot(stat_db.groupby(['ref_chr','hap'])['contig_len'].sum().reset_index(),
                x="ref_chr", y="contig_len", hue="hap")
    plt.title("len(contig)", fontsize=14)
    plt.xticks(rotation=45)
    if figName is None:
        figName = f"figs/contigLen.barplot.pdf"

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

def contigPlot(obj,width = 2, height = 4, save = True, figName = None, force = False):
    """
    Generates a heatmap of statistics for each haplotype and contig. Brick color represents T2T contigs without gaps, salmon color indicates T2T contigs with gaps, and beige color denotes non-T2T contigs.

    Parameters
    -----------
    obj
        An object that contains a `.stats` attribute, which should be a pandas DataFrame 
    """
    obj = copy.deepcopy(obj)
    stat_db = obj.stats.copy()
    stat_db.loc[stat_db['t2tStat'] == "not_t2t", "scf_ctg"] = 0
    stat_db.loc[stat_db['t2tStat'] == "scf", "scf_ctg"] = 1
    stat_db.loc[stat_db['t2tStat'] == "ctg", "scf_ctg"] = 2
    
    # Create the pivot table
    ctg = pd.pivot_table(stat_db,values='scf_ctg',index='ref_chr',columns='hap',aggfunc='max')
    # Custom labels dictionary
    custom_labels = {0: 'Not T2T', 1: 'T2T w/ gap', 2: 'T2T wo gap'}

    plt.figure(figsize=(width, height))  # Adjust the figure size as needed
    # Create a custom colormap
    cmap = ListedColormap(['#fff5f0','#fb694a','#67000d'])
    ax = sns.heatmap(ctg, cmap=cmap, linecolor="white", linewidths=0.005, cbar=False, vmin=0, vmax=2)
    handles = [mpatches.Patch(color=cmap(i), label=custom_labels[i]) for i in custom_labels]
    plt.legend(handles=handles, title='stats', loc='upper left',  bbox_to_anchor=(1, 1), frameon=False)

    # Display the plot
    if figName is None:
        figName = f"figs/contigPlot.heatmap.pdf"

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
            plt.savefig(figName, bbox_inches='tight', pad_inches=0.1    )
            print(f"File {figName} saved")

    plt.show()

    