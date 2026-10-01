import warnings
# Ensure Biopython is imported correctly
try:
    from Bio import BiopythonWarning
except ImportError:
    print("Biopython is not installed. Please install it using 'pip install biopython'.")
    sys.exit(1)

warnings.simplefilter('ignore', BiopythonWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="Bio")
warnings.filterwarnings('ignore')

from ._graphAlign import graphIdx,graphAlign,extractNodeSeq
from ._chrAssign import chrAssign,showPairwiseAlign,convertRefName,mapBetweenNodes,gfaToFasta
from ._kmer import mkMeryl,calQV
from ._asm_stats import getT2T
from ._run_rm_rDNA import rmrDNA
from ._insert_gap import insertGap
from ._modiFasta import flipContig,renameContig,sortContig,filterContigs,find_flip_candidates,map_rDNA,find_gap_in_rDNA,rDNA_gap_cleaning,screen_asm,gapCleaning
from ._detect_internal_telomere import detect_internal_telomere,runTrimming
from ._cov_unphc_phc_space import build_sparse_compression_map,lift_seqs, addPadding_to_bed, make_bandage_csv

__all__ = [
    'graphIdx',
    'extractNodeSeq',
    'graphAlign',
    'gfaToFasta',
    'chrAssign',
    'mkMeryl',
    'calQV',
    'getT2T',
    'rmrDNA',
    'insertGap',
    'convertRefName',
    'showPairwiseAlign',
    'flipContig',
    'runTrimming',
    'detect_internal_telomere',
    'renameContig',
    'sortContig',
    'filterContigs',
    'build_sparse_compression_map',
    'lift_seqs',
    'addPadding_to_bed',
    'make_bandage_csv',
    'mapBetweenNodes',
    'find_flip_candidates',
    'map_rDNA',
    'find_gap_in_rDNA',
    'rDNA_gap_cleaning',
    'screen_asm',
    'gapCleaning'
]
