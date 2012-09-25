'''
Created on Mar 9, 2012

@author: glareprotector
'''
from manager import *
import string
#import Bio.PDB
import csv
import constants
import string
import re
import math

import pdb


FILE_MANAGER_SIZE = 500
OBJ_MANAGER_SIZE = 500
MUSCLE_PATH = '/mnt/work/fultonw/muscle3.8.31_i86linux64'
DSSP_PATH = '/mnt/work/fultonw/active_site/dssp/dssp-2.0.4-linux-amd64'

to_reindex = True
recalculate = False

recalculate_nodewise_loss_f = True

metric_cutoffs = [1,2,3,4,5,6,7,8,9]




RESULTS_FOLDER = '../results/'


NACCESS_PATH = '/mnt/work/fultonw/active_site/Naccess/naccess'
NACCESS_FOLDER = '/mnt/work/fultonw/active_site/Naccess/'

BLAST_PATH = '/mnt/work/fultonw/blast/ncbi-blast-2.2.26+/bin/psiblast'
BLASTDB_PATH = '/mnt/work/fultonw/blast/ncbi-blast-2.2.26+/bin/nr/nr'
CONSERVATION_FOLDER = '/home/fultonw/conservation_code/'

ORIG_CHAINS = '../catres_pdbs'
CSA_FILE = '../catres_sites'

success_file = 'success_catres.txt'
fail_file = 'fail_catres.txt'

proc_id = 0
