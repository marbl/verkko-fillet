import pandas as pd
import os
import copy
from natsort import natsorted

def readChr(obj,mapFile, 
            chromosome_assignment_directory="chromosome_assignment", 
            stat_directory="stats",
            haplotype1=None,
            haplotype2=None):
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
        The name of haplotype 1 (e.g., "paternal" or "maternal"). Default is None.
    haplotype2
        The name of haplotype 2 (e.g., "paternal" or "maternal"). Default is None.
    
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
    translation_hap1 = pd.read_csv(f"{chromosome_assignment_directory}/translation_hap1", header = None, sep = '\t')
    translation_hap1.columns = ['contig','ref_chr','contig_len','ref_chr_len']
    translation_hap1['hap'] = translation_hap1['contig'].str.split('-', expand=True)[0]
    translation_hap1['hap'] = translation_hap1['hap'].str.split('_', expand=True)[0]
    hap1 = translation_hap1['hap'][0]
    
    translation_hap2 = pd.read_csv(f"{chromosome_assignment_directory}/translation_hap2", header = None, sep = '\t')
    translation_hap2.columns = ['contig','ref_chr','contig_len','ref_chr_len']
    translation_hap2['hap'] = translation_hap2['contig'].str.split('-', expand=True)[0]
    translation_hap2['hap'] = translation_hap2['hap'].str.split('_', expand=True)[0]
    hap2 = translation_hap2['hap'][0]
    translation = pd.concat([translation_hap1,translation_hap2])
    
    del translation_hap1
    del translation_hap2
    
    # read map filfe
    chrom_map = pd.read_csv(mapFile, sep= '\t',header = None)
    chrom_map.columns = ['old_chr','ref_chr']
 
    # read completeness 
    chr_completeness_max_hap1 = pd.read_csv(f"{chromosome_assignment_directory}/chr_completeness_max_hap1", header = None, sep = '\t')
    chr_completeness_max_hap1.columns = ['ref_chr','completeness']
    chr_completeness_max_hap1['hap']=hap1
    chr_completeness_max_hap2 = pd.read_csv(f"{chromosome_assignment_directory}/chr_completeness_max_hap2", header = None, sep = '\t')
    chr_completeness_max_hap2.columns = ['ref_chr','completeness']
    chr_completeness_max_hap2['hap']=hap2
    chr_completeness_max = pd.concat([chr_completeness_max_hap1,chr_completeness_max_hap2])
    translation['hap'] = translation['contig'].str.split('-', expand=True)[0]
    translation['hap'] = translation['hap'].str.split('_', expand=True)[0]
    
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
    
    # assembly_t2t['scf_ctg'] = pd.Categorical(assembly_t2t['scf_ctg'], categories = ["not_t2t","scf","ctg"], ordered = True)
    assembly_t2t = assembly_t2t.groupby('contig')['scf_ctg'].max()
    assembly_t2t = pd.DataFrame(assembly_t2t).reset_index()
    del assembly_t2t_scfs 
    del assembly_t2t_ctgs
    
    # merge result 
    stat_db = pd.merge(pd.merge(
        pd.merge(assembly_t2t,translation,on="contig", how = 'outer'),
        chrom_map, on='ref_chr'), chr_completeness_max, on=['ref_chr','hap'], how = 'outer')
    
    # Convert category to number
    stat_db['scf_ctg'] = stat_db['scf_ctg'].fillna(0)
    # stat_db['scf_ctg'] = pd.Categorical(stat_db['scf_ctg'], categories = ["not_t2t","scf","ctg"], ordered = True)
    
    stat_db['ref_chr'] = pd.Categorical(stat_db['ref_chr'], categories=chrom_map['ref_chr'],ordered=True)
    
    # Assuming `sire` and `dam` are defined variables
    if haplotype1!=None and haplotype2!=None:
        hap1, hap_name1 = haplotype1.split(":") 
        hap2, hap_name2 = haplotype2.split(":")
        stat_db['hap_verkko'] = stat_db['hap']
        stat_db.loc[stat_db['hap_verkko'] == hap1, "hap"] = hap_name1
        stat_db.loc[stat_db['hap_verkko'] == hap2, "hap"] = hap_name2
    else:
        stat_db['hap_verkko'] = stat_db['hap']
        print("Haplotype names not provided, using default haplotype identifiers.")
    
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