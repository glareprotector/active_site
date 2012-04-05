'''
Created on Mar 9, 2012

@author: glareprotector
'''
from manager import *
import string

FILE_MANAGER_SIZE = 500
OBJ_MANAGER_SIZE = 500




the_obj_manager = obj_manager()
the_file_manager = file_manager()

BLAST_PATH = 'psiblast'
CONSERVATION = '/home/fultonw/conservation_code'

def get_transpose(mat):
    height = len(mat)
    width = len(mat[0])
    m = [ [-1 for i in range(height)] for j in range(width)]
    for i in range(width):
        for j in range(length):
            m[i][j] = mat[j][i]
    return m

def dict_deep_copy(d):
    to_return = {}
    for key in d.keys():
        to_return[key] = d[key]
    return to_return

def write_mat(mat, f_name, the_sep = ','):
    f = open(f_name, 'w')
    for row in mat:
        line = string.join(row, sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()
