'''
Created on Mar 9, 2012

@author: glareprotector
'''
from manager import *
import string
import Bio.PDB

FILE_MANAGER_SIZE = 500
OBJ_MANAGER_SIZE = 500
MUSCLE_PATH = '/home/fultonw/muscle3.8.31_i86linux64'



the_obj_manager = obj_manager()
the_file_manager = file_manager()

BLAST_PATH = 'psiblast'
CONSERVATION_FOLDER = '/home/fultonw/conservation_code/'

ORIG_CHAINS = '../catres_pdbs'
CSA_FILE = '../catres_sites'

success_file = 'success_catres.txt'
fail_file = 'fail_catres.txt'

def get_object(p_wrapper, params, recalculate = False, to_pickle = True, use_pickle = True):
    return the_obj_manager.get_variable(p_wrapper(params), recalculate, to_pickle, use_pickle)

def get_file(p_wrapper, params, recalculate = False, option = 'r'){
    return the_file_manager.get_file(p_wrapper(params), recalculate, option)

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

def write_mat(mat, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    #print mat
    for row in mat:
        
        line = string.join([str(x) for x in row], sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()

def get_representative_atom(res):
    if 'CA' in res.child_dict.keys():
        return res['CA']
    elif 'CB' in res.child_dict.keys():
        return res['CB']
    else:
        print 'no CA or CB atom in residue'
            #pdb.set_trace()
        return res.child_list[0]
