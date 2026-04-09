import pandas as pd
import os
import re
import numpy as np
import copy


def find_intra_telo(obj, telo_file="internal_telomere/assembly_1/assembly.windows", out_prefix = None, telolen=0, loc_from_end=15000, teloPerct = 0.5):
    """\
    Find the internal telomeres in the assembly.

    Parameters
    ----------
    obj
        The object containing the stats database. 
    telo_file
        The path to the telomere file. Default is "internal_telomere/assembly_1/assembly.windows". 
    out_prefix
        The prefix for the output file. Default is None. 
    telolen
        The minimum length of the telomere. Default is 0. 
    loc_from_end
        The minimum distance from the end of the contig. Default is 15000.  
    teloPerct
        The minimum telomere percentage. Default is 0.5. 

    Returns
    -------
        The DataFrame containing the contig, internal-left, internal-right, distal-left, distal-right, and problem columns.
    """
    print("Finding the internal telomeres in the assembly...")

    obj = copy.deepcopy(obj)
    # Check if the stats database is available
    if not isinstance(obj.stats, pd.DataFrame):
        raise ValueError("The stats database is not available. Please run the readChr function first.")

    statsdb = obj.stats

    telo_file = os.path.abspath(telo_file)
    if not os.path.exists(telo_file):
        raise FileNotFoundError(f"File not found: {telo_file}")
    print(f"Reading file: {telo_file}")
    tel = pd.read_csv(telo_file, sep='\t', header=None, usecols=[1, 2, 3, 4, 5])
    tel.columns = ['contig', 'totalLen', 'start', 'end', 'teloPerct']
    tel['contig'] = tel['contig'].str.replace('^>', '', regex=True)
    tel['start'] = tel['start'].astype(int)
    tel['end'] = tel['end'].astype(int)
    tel['totalLen'] = tel['totalLen'].astype(int)
    # Filter based on teloPerct
    tel = tel[tel['teloPerct'] > teloPerct]
    tel.reset_index(drop=True, inplace=True)

    # Iterate over rows to merge them
    rows_to_remove = []  # To keep track of rows to be dropped after merging
    for i in range(len(tel) - 1, 0, -1):  # Iterate backward to avoid index issues
        if (tel.loc[i, 'start'] < tel.loc[i - 1, 'end']) and \
        (tel.loc[i, 'contig'] == tel.loc[i - 1, 'contig']):
            
            # Merge the rows by combining the 'start', 'end', and other columns
            tel.loc[i - 1, 'end'] = tel.loc[i, 'end']  # Update 'end' to the next row's 'end'
            tel.loc[i - 1, 'teloPerct'] = max(tel.loc[i, 'teloPerct'], tel.loc[i - 1, 'teloPerct'])
            # Mark the current row for removal (row i)
            rows_to_remove.append(i)
    # Drop the rows that have been merged
    tel = tel.drop(rows_to_remove).reset_index(drop=True)


    tel['telomere'] = 'distal'
    tel['contig_start_len'] = tel['start']
    tel['contig_end_len'] = tel['totalLen'] - tel['end']
    tel['lenmin'] = tel[['contig_start_len', 'contig_end_len']].min(axis=1)
    tel.loc[tel['lenmin'] > loc_from_end, 'telomere'] = 'internal'
    intTelCount = tel[tel['telomere'] == 'internal'].shape[0]
    print(f"Number of internal telomeres: {intTelCount}")

    tel['arm'] = ""
    tel.loc[tel['start'] < tel['totalLen'] / 2, 'arm'] = "left"
    tel.loc[tel['start'] > tel['totalLen'] / 2, 'arm'] = "right"

    # Filter based on length (if needed)
    tel['len'] = tel['end'] - tel['start']
    tel = tel.loc[tel['len'] >= telolen]
    
    # Filter contigs that are not in the stats database
    tel = tel[tel['contig'].isin(obj.stats['contig'])]
    # drop coullmns
    tel.drop(columns = ['lenmin', 'contig_start_len','contig_end_len'], inplace=True)

    tel['tel-arm'] = tel['telomere'] + "-" + tel['arm']
    result = tel.groupby(['contig','tel-arm'])['teloPerct'].max().unstack(fill_value=0)
    result['problem'] = "OK/OK"
    # check if columns are in the result
    check_columns = ['internal-left', 'internal-right', 'distal-left', 'distal-right']
    for col in check_columns:
        if col not in result.columns:
            result[col] = 0

    missingTel = (result['distal-left']==0) | (result['distal-right']==0) 
    INTEL = (result['internal-left']==0) | (result['internal-right']==0)

    result.loc[missingTel & INTEL, 'problem'] = "MissingTel/INTEL"
    result.loc[missingTel & ~INTEL, 'problem'] = "MissingTel/OK"
    result.loc[~missingTel & INTEL, 'problem'] = "OK/INTEL"

    result_merged = pd.merge(result, statsdb, left_on='contig', right_on='contig', how='right')
    
    
    # result_merged
    if out_prefix is None:
        out_prefix = f"{telo_file}.loc_from_end_{loc_from_end}.teloPerct_{teloPerct}.telolen_{telolen}.stats"
    wirte_to = out_prefix + ".tsv"
    result_merged.to_csv(wirte_to, sep='\t', index=False)
    print(f"File saved: {wirte_to}")
    
    return result_merged, tel


def find_reads_intra_telo(tel, lineNum ,scfmap = "assembly.scfmap",layout = "6-layoutContigs/unitig-popped.layout", gap_bed= "stats/assembly.gaps.bed"):
    """\
    Find the reads support for the additional artifical sequences outside of the telomere.

    Parameters
    ----------
    tel
        The DataFrame containing the telomere information. This is the output of the find_intra_telo function.
    lineNum
        The line number of the telomere. This is the index of the tel DataFrame.
    scfmap
        The path to the scfmap file. Default is "assembly.scfmap". 
    layout
        The path to the layout file. Default is "6-layoutContigs/unitig-popped.layout". 

    Returns
    -------
        The DataFrame containing the readName, start_hpc, end_hpc, start, end, and type columns.
    """
    print("Finding the reads support for the additional artifical sequences outside of the telomere...")

    intra_telo = tel.copy()    
    intra_telo = intra_telo.loc[lineNum,:]
    contig = intra_telo['contig']

    if (intra_telo['start'] - 0) > (intra_telo['totalLen']-intra_telo['end']):
        pos= "end"
    else:
        pos = "start"
    print(f"Looking for the reads from {pos} of {contig}")

    if pos == 'start':
        bp = int(intra_telo['start'])
    elif pos == 'end' :
        bp = int(intra_telo['end'])
    else :
        print ("the pos argument should be either start or end")
        return
    len_fai = int(intra_telo['totalLen'])

    with open(scfmap, 'rb') as f:
        data = f.read().decode('utf-8')  # Decode bytes to string
    
    # Regular expression to match 'path' to 'end'
    pattern = r'(path.*?end)'
    
    # Find all matches
    matches = re.findall(pattern, data, re.DOTALL)
    filtered_matches = [match for match in matches if contig in match]
    
    match_pieces = filtered_matches[0].split("\n") if filtered_matches else []
    match_pieces.pop(0)
    match_pieces.pop(-1)
    # match_pieces

    # Get the first and last piece
    if pos == 'start':
        piece_selected = match_pieces[0] if match_pieces else None
    elif pos == 'end':
        piece_selected = match_pieces[-1] if match_pieces else None

    if len(match_pieces) > 1 : 
        gap = pd.read_csv(gap_bed, sep="\t", header=None, names=["chrom", "start", "end"])
        gap_sub = gap.loc[gap['chrom'] == contig, ]
        if pos == 'start':
            piece_start = 0
            piece_end = gap_sub.loc[gap_sub.index[0], 'start'] -1
        if pos == 'end':
            piece_start = gap_sub.loc[gap_sub.index[-1], 'end'] + 1
            piece_end = len_fai

    else:
        piece_start = 0
        piece_end = len_fai
    print(f"the piece of interest is from {piece_start} to {piece_end}, match piece is {piece_selected}")

    with open(layout, 'rb') as f:
        data = f.read().decode('utf-8')  # Decode bytes to string

    # Regular expression to match 'path' to 'end'
    pattern = r'(tig.*?end)'

    # Find all matches
    matches = re.findall(pattern, data, re.DOTALL)
    filtered_matches = [match for match in matches if piece_selected in match]
    filtered_matches = filtered_matches[0].split("\n")

    filtered_matches_body = filtered_matches[4:-1]
    filtered_matches_body = [entry.split("\t") for entry in filtered_matches_body]
    # print(filtered_matches_body)
    df = pd.DataFrame(filtered_matches_body, columns=["readName", "start_hpc", "end_hpc","nTime","Quality"])
    df = df[['readName', 'start_hpc', 'end_hpc','Quality']]

    df['end_hpc'] = df['end_hpc'].astype(int)
    df['start_hpc'] = df['start_hpc'].astype(int)
    piece_len_hpc = max(df['start_hpc'].max(), df['end_hpc'].max())

    piece_len = piece_end - piece_start

    print(f"the piece length is {piece_len}, the piece length in hpc is {piece_len_hpc}")

    df['start_hpc'] = df['start_hpc'].astype(int)
    df['end_hpc'] = df['end_hpc'].astype(int)

    # scaling the hpc to real length
    df['start'] = (df[['end_hpc', 'start_hpc']].min(axis=1) * piece_len / piece_len_hpc) + piece_start
    df['end'] = (df[['end_hpc', 'start_hpc']].max(axis=1) * piece_len / piece_len_hpc ) + piece_start
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

    pieceinfo = filtered_matches[0:4]
    pieceinfo = [entry.split("\t") for entry in pieceinfo]
    
    def classify(name):
        if name.startswith('ont'):
            return 'ont'
        elif name.startswith('hifi_m'):
            return 'hifi'
        elif name.startswith('hifi'):
            return 'ont_corrected'
        else:
            return None

    df['type'] = df['readName'].apply(classify)
    df['type'] = pd.Categorical(df['type'], categories=['ont','hifi','ont_corrected'], ordered=True)

    if pos == "start":
        # df_sub = df.loc[(df['start'] < bp)|(df['end'] < bp)]
        df_sub = df.loc[df['start'].between(0,bp)]
    elif pos == "end":
        # df_sub = df.loc[(df['start'] > bp)|(df['end'] > bp)]
        df_sub = df.loc[df['end'].between(bp,len_fai)]
    else:
        print("pos should be either start or end")
    
    df_sub_count = df_sub.groupby('type')['start_hpc'].count().reset_index()
    
    print("Summary : ")
    print("   Num of ONT reads : " + str(df_sub_count.iloc[0,1]))
    print("   Num of HiFi reads : " + str(df_sub_count.iloc[1,1]))
    print("   Num of ONT corrected reads : " + str(df_sub_count.iloc[2,1]))
    
    return df_sub, df