import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from .._run_shell import run_shell
import os
import subprocess
import matplotlib
from natsort import natsorted
import copy
from PIL import Image
from IPython.display import display
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def handle_duplicates(series):
    counts = {}  # To keep track of occurrences of each value
    new_values = []
    for value in series:
        if value in counts:
            counts[value] += 1
            new_values.append(f"{value}_{counts[value]}")
        else:
            counts[value] = 1
            new_values.append(value)
    return new_values


def percTel(intra_telo , showContig=None,
            width = 5, height = 7, save = True, figName = None):
    """
    Generates a heatmap showing the telomere percentage by contig.

    Parameters
    -----------
    intra_telo
        A DataFrame containing the telomere percentage data. 
    showContig
        Columns to show in the heatmap. Default is None. If None, only the 'contig' column is shown. 
    width
        Width of the plot. Default is 5. 
    height
        Height of the plot. Default is 7.
    save
        If True, the plot is saved as a PNG file. Default is True. 
    figName
        Name of the saved plot. Default is None. If None, the plot is saved as "figs/intra_telo.heatmap.png". 
    """
    
    if showContig is None:
        showContig = ['contig']

    check_columns = ['distal-left', 'internal-left', 'internal-right', 'distal-right']
    intra_telo = intra_telo.copy()
    intra_telo['by'] = intra_telo[showContig].astype(str).agg('_'.join, axis=1)
    heatmapDb = intra_telo.loc[:,['by']+ check_columns]
    heatmapDb['by'] = handle_duplicates(heatmapDb['by'])
    heatmapDb = heatmapDb.set_index('by')
    heatmapDb = heatmapDb.dropna()
    heatmapDb_idx = natsorted(heatmapDb.index)
    heatmapDb = heatmapDb.loc[heatmapDb_idx, check_columns]


    plt.figure(figsize=(width, height))

    ax = sns.heatmap(heatmapDb, annot=True, fmt=".2f", cmap = "Reds", 
                    cbar_kws={'label': 'Telomere Percentage'}, vmin=0, vmax=1,  linewidth=.3)
    ax.set(xlabel="", ylabel="")
    ax.set(title="Telomere Percentage by Contig")

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




def readOnNode(tel, lineNum, readBed, figName = None):
    tel = tel.copy()
    df_copy = readBed.copy()
    ind = lineNum

    # Ensure directory exists
    os.makedirs("internal_telomere", exist_ok=True)


    df_copy['start'] = df_copy['start'].astype(int)
    df_copy['end'] = df_copy['end'].astype(int)
    df_copy['type'] = df_copy['type'].astype(str)
    df_copy['start'] = df_copy['start'].astype(int)
    df_copy['end'] = df_copy['end'].astype(int)

    # tel = tel.copy()  # Ensure tel is defined earlier
    contig = tel.loc[ind, "contig"].replace(".", "_").replace("-", "_")
    telomere  = tel.loc[ind, "telomere"]
    arm = tel.loc[2, "arm"]
    prefix = f"{contig}_{telomere}_{arm}"
    print(f"Prefix: {prefix}")


    if arm == "left":
        show_start = 1
        show_end = tel.loc[ind, "end"]
    else:
        show_start = tel.loc[ind, "start"]
        show_end = tel.loc[ind, "totalLen"]
    print(f"show region: {contig}:{show_start}-{show_end}")

    # vbedf 
    vhighlightBed = f"internal_telomere/{prefix}.v.bed"
    if os.path.exists(vhighlightBed):
        print("Highlight BED file exists")
    else:
        print("Creating highlight BED file")
        v_start= tel.loc[ind, "start"]
        v_end = tel.loc[ind, "end"]
        vbed = pd.DataFrame({"chrom": [contig], "start": [v_start], "end": [v_end]})
        vbed.to_csv(vhighlightBed, sep="\t", header=False, index=False)
        print(f"highlight region: {contig}:{v_start}-{v_end}")

    # Prepare BED file from df (ensure df has chrom, start, end columns)
    readBed = f"internal_telomere/{prefix}.bed"
    if os.path.exists(readBed):
        print("Read BED file exists")
    else:
        
        print("Creating BED file")
        df_bed = df_copy.loc[:, ['readName','start','end','type']]  # Ensure columns are correct
        df_bed.columns = ["chrom", "start", "end", 'type']
        df_bed["chrom"] = contig
        df_bed["strand"] = "+"

        df_bed[['start', 'end', 'strand']] = df_bed.apply(
            lambda x: (x['start'], x['end'], '+') if x['start'] < x['end'] else (x['end'], x['start'], '-'),
            axis=1, result_type="expand"
        )
        df_bed# Ensure columns are correct
        df_bed = df_bed.loc[:, ['chrom', 'start', 'end', 'type']]
        df_bed.columns = ["chrom", "start", "end", 'type']
        df_bed["chrom"] = contig  # Make sure 'contig' is defined
        df_bed["start"] = df_bed["start"].astype(int)
        df_bed["end"] = df_bed["end"].astype(int)
        df_bed["type"] = df_bed["type"].astype(str)

        # Apply swap and strand assignment in one step
        df_bed[['start', 'end', 'strand']] = df_bed.apply(
            lambda x: (x['start'], x['end'], '+') if x['start'] < x['end'] else (x['end'], x['start'], '-'),
            axis=1, result_type="expand"
        )
        df_bed['name'] = "0"
        df_bed['score'] = "0"
        df_bed['thickStart'] = df_bed['start']
        df_bed['thickEnd'] = df_bed['end']

        ### Color by type
        # df_bed['color'] = df_bed['type'].apply(lambda x: '255,0,0' if x == 'ont' else '128,255,0')
        def classify(name):
            if name == 'ont':
                return '255,0,0' #red
            elif name == "hifi":
                return '128,255,0' # green
            elif name == "ont_corrected":
                return "245,166,39" # orange
            else:
                return None

        df_bed['color'] = df_bed['type'].apply(classify)
        df_bed = df_bed[['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'color']]

        # Write BED file

        comment=f"track name=\"ItemRGBDemo\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\""
        subprocess.run(f"echo {comment} > {readBed}", shell=True)
        df_bed.to_csv(readBed, sep="\t", header=False, index=False, mode = 'a')

    # Track file
    trackFile = f"internal_telomere/{prefix}.tracks.ini"
    if os.path.exists(trackFile):
        print("Track file exists")
    else:
        # result = subprocess.run(
            #f"make_tracks_file --trackFiles internal_telomere/{prefix}.read.bed -o internal_telomere/{prefix}.read.bed.tracks.ini",
            #check=True, shell=True, capture_output=True, text=True)

        track_contents = [
            "[x-axis]",
            " ",
            "[spacer]",
            "height = .3",
            " ",
            f"[{contig}]",
            f"file = {readBed}",
            f"title = Composition of reads",
            "height = 30",
            "color = bed_rgb",
            "border_color = None",
            "labels = false",
            "fontsize = 10",
            "file_type = bed",
            " ",
            f"[{contig}_highlight]",
            f"file = {vhighlightBed}",
            "type = vhighlight",
            "alpha = 0.2",
            "zorder = -100"
        ]

        if os.path.exists(trackFile):
            os.remove(trackFile)
        for elements in track_contents:
            subprocess.run(f"echo {elements} >> {trackFile}", shell=True)

    
    figName=f"internal_telomere/{prefix}_output.png" if figName is None else figName

    if os.path.exists(figName):
        print("Plot file exists")
    else:
        print("Creating plot file")
        result = subprocess.run(
            f"pyGenomeTracks --tracks {trackFile} --dpi 400 --region {contig}:{show_start}-{show_end} --outFileName {figName}",
            check=True, shell=True, capture_output=True, text=True
        )
        print("pyGenomeTracks STDOUT:", result.stdout)
        print("pyGenomeTracks STDERR:", result.stderr)

    try:
        img = Image.open(figName)
    except FileNotFoundError:
        print(f"Error: Image file not found at '{figName}'")
    else:
        # Display the image
        display(img)