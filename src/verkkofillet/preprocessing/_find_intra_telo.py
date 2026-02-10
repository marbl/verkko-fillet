import pandas as pd
import os
import re
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


def find_reads_intra_telo(tel, lineNum ,scfmap = "assembly.scfmap",layout = "6-layoutContigs/unitig-popped.layout"):
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
    
    # Regular expression to match all pieces
    pattern = r'piece\d{6}'
    
    # Extract pieces from the data string
    pieces = re.findall(pattern, filtered_matches[0])
    
    # Get the first and last piece
    first_piece = pieces[0] if pieces else None
    last_piece = pieces[-1] if pieces else None
    
    # Output the results
    # print(f"First piece: {first_piece}")
    # print(f"Last piece: {last_piece}")
    
    if pos == "start":
        piece = first_piece
    elif pos == "end":
        piece = last_piece
    else:
        print("pos should be either start or end")

    print("Looking for the reads from " + piece)
    
    with open(layout, 'rb') as f:
        data = f.read().decode('utf-8')  # Decode bytes to string
    
    # Regular expression to match 'path' to 'end'
    pattern = r'(tig.*?end)'
    
    # Find all matches
    matches = re.findall(pattern, data, re.DOTALL)
    filtered_matches = [match for match in matches if piece in match]
    filtered_matches = filtered_matches[0].split("\n")

    filtered_matches_body = filtered_matches[4:-1]
    filtered_matches_body = [entry.split("\t") for entry in filtered_matches_body]
    # print(filtered_matches_body)
    df = pd.DataFrame(filtered_matches_body, columns=["readName", "start_hpc", "end_hpc","nTime","Quality"])
    df = df[['readName', 'start_hpc', 'end_hpc','Quality']]
    
    df['start_hpc'] = df['start_hpc'].astype(int)
    df['end_hpc'] = df['end_hpc'].astype(int)
    
    df['start'] = df[['end_hpc', 'start_hpc']].min(axis=1) * 1.6  # multiply 1.6 cuz this is based on HPC coordinates
    df['end'] = df[['end_hpc', 'start_hpc']].max(axis=1) * 1.6    # multiply 1.6 cuz this is based on HPC coordinates
    
    pieceinfo = filtered_matches[0:4]
    pieceinfo = [entry.split("\t") for entry in pieceinfo]

    df['type'] = df['readName'].apply(lambda x: 'ont' if 'ont' in x else 'hifi')
    df['type'] = pd.Categorical(df['type'], categories=['ont','hifi'], ordered=True)

    if pos == "start":
        df_sub = df.loc[(df['start'] < bp)|(df['end'] < bp)]
    elif pos == "end":
        df['start'] = len_fai - int(pieceinfo[1][1]) + df['start']
        df['end'] = len_fai - int(pieceinfo[1][1]) + df['end']
        df_sub = df.loc[(df['start'] > bp)|(df['end'] > bp)]
    else:
        print("pos should be either start or end")
    
    df_sub_count = df_sub.groupby('type')['start_hpc'].count().reset_index()
    
    print("Summary : ")
    print("   Num of ONT reads : " + str(df_sub_count.iloc[0,1]))
    print("   Num of HiFi reads : " + str(df_sub_count.iloc[1,1]))
    
    return df_sub, df