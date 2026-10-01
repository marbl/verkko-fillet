import pandas as pd
import os
import copy
from natsort import natsorted

def readChr(obj,mapFile, 
            chromosome_assignment_directory="chromosome_assignment", 
            stat_directory="stats",
            haplotype1="haplotype1",
            haplotype2="haplotype2"):
    """\
    Read the chromosome assignment results and store them in the object.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    mapFile
        The path to the map file.
    chromosome_assignment_directory
        The directory containing the chromosome assignment results. Default is "chromosome_assignment".
    stat_directory
        The directory containing the statistics. Default is "stats".
    haplotype1
        The name of haplotype 1 (e.g., "paternal" or "maternal"). Default is haplotype1.
    haplotype2
        The name of haplotype 2 (e.g., "paternal" or "maternal"). Default is haplotype2.
    
    Returns
    -------
        obj with stats attribute containing the chromosome assignment results.
    """
    obj = copy.deepcopy(obj)
    chromosome_assignment_directory = os.path.abspath(chromosome_assignment_directory)
    stat_directory = os.path.abspath(stat_directory)
    mapFile = os.path.abspath(mapFile)

    if not os.path.exists(chromosome_assignment_directory):
        raise FileNotFoundError(f"{chromosome_assignment_directory} does not exist.")
    if not os.path.exists(stat_directory):
        raise FileNotFoundError(f"{stat_directory} does not exist.")
    if not os.path.exists(mapFile):
        raise FileNotFoundError(f"{mapFile} does not exist.")
    
    # read translation
    print(f"Reading translation files from {chromosome_assignment_directory}")
    translation_hap1 = pd.read_csv(f"{chromosome_assignment_directory}/translation_hap1", header = None, sep = '\t')
    translation_hap1.columns = ['contig','ref_chr','contig_len','ref_chr_len']
    # check if haplotype1 is in the contig names
    if len(translation_hap1['contig'].str.contains(haplotype1)) == 0:
        raise ValueError(f"Haplotype1 {haplotype1} not found in contig names.")
    translation_hap1['hap'] = haplotype1 if translation_hap1['contig'].str.contains(haplotype1).all() else haplotype2 if translation_hap1['contig'].str.contains(haplotype2).all() else None
        
    translation_hap2 = pd.read_csv(f"{chromosome_assignment_directory}/translation_hap2", header = None, sep = '\t')
    translation_hap2.columns = ['contig','ref_chr','contig_len','ref_chr_len']
    if len(translation_hap2['contig'].str.contains(haplotype2)) == 0:
        raise ValueError(f"Haplotype2 {haplotype2} not found in contig names.")
    translation_hap2['hap'] = haplotype2 if translation_hap2['contig'].str.contains(haplotype2).all() else haplotype1 if translation_hap2['contig'].str.contains(haplotype1).all() else None
    
    translation = pd.concat([translation_hap1,translation_hap2])
    
    del translation_hap1
    del translation_hap2
    
    # read map filfe
    print(f"Reading map file from {mapFile}")
    chrom_map = pd.read_csv(mapFile, sep= '\t',header = None)
    chrom_map.columns = ['old_chr','ref_chr']
 
    # read completeness 
    print(f"Reading completeness files from {chromosome_assignment_directory}")
    chr_completeness_max_hap1 = pd.read_csv(f"{chromosome_assignment_directory}/chr_completeness_max_hap1", header = None, sep = '\t')
    chr_completeness_max_hap1.columns = ['ref_chr','completeness']
    chr_completeness_max_hap1['hap'] = haplotype1
    chr_completeness_max_hap2 = pd.read_csv(f"{chromosome_assignment_directory}/chr_completeness_max_hap2", header = None, sep = '\t')
    chr_completeness_max_hap2.columns = ['ref_chr','completeness']
    chr_completeness_max_hap2['hap']=haplotype2

    chr_completeness_max = pd.concat([chr_completeness_max_hap1,chr_completeness_max_hap2])
    del chr_completeness_max_hap2
    del chr_completeness_max_hap1
    
    # read t2t stat
    assembly_t2t_scfs = pd.read_csv(f"{stat_directory}/assembly.t2t_scfs", header = None, sep = '\t')
    assembly_t2t_scfs.columns = ['contig']
    assembly_t2t_scfs['scf_ctg'] = 1
    assembly_t2t_ctgs = pd.read_csv(f"{stat_directory}/assembly.t2t_ctgs", header = None, sep = '\t')
    assembly_t2t_ctgs.columns = ['contig']
    assembly_t2t_ctgs['scf_ctg'] = 2
    assembly_t2t = pd.concat([assembly_t2t_scfs,assembly_t2t_ctgs])
    assembly_t2t = assembly_t2t.groupby('contig')['scf_ctg'].max()
    assembly_t2t = pd.DataFrame(assembly_t2t).reset_index()
    del assembly_t2t_scfs 
    del assembly_t2t_ctgs
    
    # merge result 
    print("Merging results...")
    stat_db = pd.merge(pd.merge(
        pd.merge(assembly_t2t,translation,on="contig", how = 'outer'),
        chrom_map, on='ref_chr'), chr_completeness_max, on=['ref_chr','hap'], how = 'outer')

    # Convert category to number
    stat_db['scf_ctg'] = stat_db['scf_ctg'].fillna(0)
    stat_db['ref_chr'] = pd.Categorical(stat_db['ref_chr'], categories=chrom_map['ref_chr'],ordered=True)
    
    stat_db.loc[stat_db['scf_ctg'] == 0, "t2tStat"] = "not_t2t"
    stat_db.loc[stat_db['scf_ctg'] == 1, "t2tStat"] = "scf"
    stat_db.loc[stat_db['scf_ctg'] == 2, "t2tStat"] = "ctg"
    del stat_db['scf_ctg']
    
    obj.stats = stat_db
    print("The chromosome infomation was stored in obj.stats")
    return obj

def detectBrokenContigs(obj):
    """\
    Find contigs that assigned same chromosome and haplotype.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.

    Returns
    -------
        The DataFrame containing the duplicated contigs with different contig names.
    """
    obj = copy.deepcopy(obj)
    stats = obj.stats.copy()
    tab = stats.groupby(["ref_chr","hap"])['contig'].count().reset_index()
    tab = tab.loc[tab['contig']>1].reset_index(drop=True)
    if len(tab)>0:
        print("Warning: the following chromosomes have more than one contig:")
        print(tab)
    else:
        print("All chromosomes have one contig!")