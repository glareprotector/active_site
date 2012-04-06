'''
Created on Mar 9, 2012

@author: glareprotector
'''
from manager import *
import string

FILE_MANAGER_SIZE = 500
OBJ_MANAGER_SIZE = 500
MUSCLE_PATH = '/home/fultonw/muscle3.8.31_i86linux64'



the_obj_manager = obj_manager()
the_file_manager = file_manager()

BLAST_PATH = 'psiblast'
CONSERVATION_FOLDER = '/home/fultonw/conservation_code/'

ORIG_CHAINS = '../chains_to_use.txt'
CSA_FILE = '../CSA_2_2_12.dat'

def get_transpose(mat):
    height = len(mat)
    width = len(mat[0])
    m = [ [-1 for i in range(height)] for j in range(width)]
    for i in range(width):
        for j in range(height):
            m[i][j] = mat[j][i]
    return m

def dict_deep_copy(d):
    to_return = {}
    for key in d.keys():
        to_return[key] = d[key]
    return to_return

def write_mat(mat, f_name, the_sep = ','):
    f = open(f_name, 'w')
    #print mat
    for row in mat:
        
        line = string.join([str(x) for x in row], sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()
