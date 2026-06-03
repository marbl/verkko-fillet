#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys 
import pandas as pd
import os
pd.set_option('mode.chained_assignment', None)
import warnings
import session_info
warnings.simplefilter(action='ignore', category=FutureWarning)
import verkkofillet as vf


# In[ ]:
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--ref",        help="reference fasta file")
parser.add_argument("--map_file",   help=(
    "chrMap file: tab-separated mapping from contig names in the reference\n"
    "fasta to the desired chromosome names.\n"
    "  Format : <contig_name_in_fasta>\t<desired_chromosome_name>\n"
    "  Examples:\n"
    "    chr1\tchr1\t# contig is already named as desired\n"
    "    GCF00001\tchr1\t# rename contig GCF00001 to chr1"
))
parser.add_argument("--internal_tel_threshold", type=float, help="internal telomere threshold, default is 15000 bp", default=15000)
parser.add_argument("--verkkoDir",  help="path to verkko output directory")
args, unknown = parser.parse_known_args()

ref = args.ref
map_file = args.map_file
internal_tel_threshold = args.internal_tel_threshold
verkkoDir = args.verkkoDir

missing = [name for name, value in (("--verkkoDir", verkkoDir), ("--ref", ref), ("--map_file", map_file)) if not value]
if missing:
    print(f"Error: {', '.join(missing)} argument(s) are required.")
    sys.exit(1)

if not os.path.exists(ref):
    print(f"Error: Reference fasta file '{ref}' does not exist.")
    sys.exit(1)
if not os.path.exists(map_file):
    print(f"Error: Map file '{map_file}' does not exist.")
    sys.exit(1)
if not os.path.exists(verkkoDir):
    print(f"Error: Verkko directory '{verkkoDir}' does not exist.")
    sys.exit(1)


print(f"Reference: {ref}")
print(f"Map file: {map_file}")
print(f"Internal telomere threshold: {internal_tel_threshold}")
print(f"Verkko directory: {verkkoDir}")

if not os.path.exists(verkkoDir):
    print(f"Error: Verkko directory '{verkkoDir}' does not exist.")
    sys.exit(1)

os.chdir(verkkoDir)

# In[3]:


obj = vf.pp.read_Verkko(verkkoDir, lock_original_folder = False)
os.getcwd()

# In[4]:


obj


# In[13]:



# In[17]:

vf.tl.getT2T(obj)


vf.tl.chrAssign(obj = obj, ref = ref)


# In[5]:


obj = vf.pp.readChr(obj, map_file)
obj.stats


# In[6]:


vf.pp.detectBrokenContigs(obj)


# In[7]:


vf.pl.showMashmapOri(obj, height = 5 , width = 5)


# In[8]:


vf.pl.completePlot(obj, height = 3 , width = 6)


# In[9]:


vf.pl.contigLenPlot(obj,height = 3 , width = 6)


# In[10]:


vf.pl.contigPlot(obj,height = 6 , width = 2)


# In[11]:


vf.pl.n50Plot(obj, height = 4 , width = 10)


# In[ ]:


vf.tl.detect_internal_telomere(obj)
result_merged, tel = vf.pp.find_intra_telo(obj, loc_from_end= internal_tel_threshold)
result_merged


# In[ ]:


vf.pl.percTel(result_merged, showContig  = ['ref_chr','hap'], height = 10 , width = 5)

print("Done!")
print(f"All figures are saved in {obj.verkko_fillet_dir}/figs/ folder.")
print(f"All stats are saved in {obj.verkko_fillet_dir}/stats/ folder.")
print(f"All chromosome assignment files are saved in {obj.verkko_fillet_dir}/chromosome_assignment/ folder.")

# session info
print("\nSession Info:")
session_info.show()