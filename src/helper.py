import string
#import Bio.PDB
import csv
import constants
import string
import re
import math
import global_stuff
import pdb
import param

# also sets global_stuff.RESULTS_FOLDER to proper value
def read_param(file_location):
    
    # read folder_name
    f = open(constants.INFO_FOLDER + file_location)
    the_params = param.param({})
    hp_values = param.param()
    folder_name = f.readline().strip()

    global_stuff.RESULTS_FOLDER = global_stuff.RESULTS_BASE_FOLDER + folder_name + '/'

    for line in f.readlines():
        print line
        if line[0] != '#':
            s = line.strip().split(',')
            if s[0] != 'hp':
                the_name = s[0]
                if the_name == 'n':
                    node_features = []
                    for i in range(1, len(s)):
                        node_features.append(constants.get_master_node_feature_list()[int(s[i])])
                    the_params.set_param('n', node_features)
                if the_name == 'e':
                    edge_features = []
                    for i in range(1, len(s)):
                        edge_features.append(constants.get_master_edge_feature_list()[int(s[i])])
                    the_params.set_param('e', edge_features)
                try:
                    the_type = s[1]
                    if the_type == 'f':
                        the_params.set_param(the_name, float(s[2]))
                    elif the_type == 'i':
                        the_params.set_param(the_name, int(s[2]))
                    elif the_type == 's':
                        the_params.set_param(the_name, s[2])
                except:
                    pass

    # hp values file happens to be the same as info file, so set that

    the_params.set_param('hpvf', file_location)


    if len(the_params.get_param('e')) != 0:
        assert the_params.get_param('wif') != 2

    
    return folder_name, the_params

def read_hp_values(file_location):

    # read folder_name
    f = open(constants.INFO_FOLDER + file_location)

    hp_values = param.param()
    folder_name = f.readline().strip()



    for line in f.readlines():
        s = line.strip().split(',')
        if s[0] == 'hp':
            to_add = []
            the_type = s[2]
            the_name = s[1]
            for i in range(3,len(s)):
                if the_type == 'f':
                    to_add.append(float(s[i]))
                elif the_type == 'i':
                    to_add.append(int(s[i]))
                elif the_type == 's':
                    to_add.append(s[i])
            hp_values.set_param(the_name, to_add)

    return hp_values

def get_aux_folder(pdb_name, chain_letter, start, end):
    return constants.AUX_FOLDER + string.lower(pdb_name) + '_' + string.upper(chain_letter) + '_' + str(start) + '_' + str(end) + '/'

# returns mat normalized by columns
def normalize_mat(mat):
    height = len(mat)
    width = len(mat[0])
    normalized = [ [0 for i in range(width)] for j in range(height)]

    for i in range(width):
        a_mean_sum = 0.0
        for j in range(height):
            a_mean_sum = a_mean_sum + mat[j][i]
        a_mean = a_mean_sum / height
        a_var_sum = 0
        for j in range(height):
            a_var_sum = a_var_sum + pow((mat[j][i] - a_mean), 2)
        sd = math.sqrt(a_var_sum / height)
        if abs(sd) > .00001:
            for j in range(height):
                normalized[j][i] = (mat[j][i] - a_mean) / sd
    return normalized


def shorten(x):
    x = re.sub(r'\'','',x)
    x = re.sub(r'class','',x)
    x = re.sub(r' ','',x)
    x = re.sub(r'<','',x)
    x = re.sub(r'>','',x)
    x = re.sub(r'f\.','',x)
    x = re.sub(r'\),\(',')(',x)
    return x


def super_shorten(x):

    x = re.sub(r'\'','',x)
    x = re.sub(r'class','',x)
    x = re.sub(r' ','',x)

    
    x = re.sub(r'<','',x)
    x = re.sub(r'>','',x)
    x = re.sub(r'f\.','',x)
    x = re.sub(r'\),\(',')(',x)
    #x = re.sub(r'\)\(','|',x)
    #x = re.sub(r'\[\(','[',x)
    #x = re.sub(r'\)\]',']',x)
    #pdb.set_trace()
    return x


def get_KL(d1, d2):

    # keep dictionary of counts for each distribution
    d1_dict = {}
    for i in range(len(d1)):
        if d1[i] in d1_dict.keys():
            d1_dict[d1[i]] = d1_dict[d1[i]] + 1
        else:
            d1_dict[d1[i]] = 1
    for k in d1_dict.keys():
        d1_dict[k] = float(d1_dict[k]) / float(len(d1))
 #   pdb.set_trace()
    d2_dict = {}
    for i in range(len(d2)):
        if d2[i] in d2_dict.keys():
            d2_dict[d2[i]] = d2_dict[d2[i]] + 1
        else:
            d2_dict[d2[i]] = 1
    for k in d2_dict.keys():
        d2_dict[k] = float(d2_dict[k]) / float(len(d2))
#    pdb.set_trace()
    # add pseudocounts for d2
    for k in d1_dict.keys():
        if k in d2_dict.keys():
            d2_dict[k] = d2_dict[k] + 1
        else:
            d2_dict[k] = 1
    ans = 0
    for k in d1_dict.keys():
        if d1_dict[k] > 0:
            ans += d1_dict[k] * math.log(d1_dict[k] / d2_dict[k])
    return ans


def physical_distance(x):
    norm = 0;
    for i in range(len(x)):
        norm = norm + x[i] * x[i]
    return math.sqrt(norm)

def print_stuff_dec(f):

    def g(*args, **kwargs):
        #print 'calling ', f.func_name, ' with ', args, kwargs
        ans = f(*args, **kwargs)
        #print f.func_name, ' returned ', ans
        return ans
    
    return g

def get_object(p_wrapper, params, recalculate = False, to_pickle = True, use_pickle = True):
    return the_obj_manager.get_variable(p_wrapper(params), recalculate, to_pickle, use_pickle)

def get_file(p_wrapper, params, recalculate = False, option = 'r'):
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

def list_union(a, b):
    A = set(a)
    B = set(b)
    return list(A - B)


def write_mat(mat, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    #print mat
    for row in mat:
        
        line = string.join([str(x) for x in row], sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()

def write_vect(vect, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    line = string.join([str(x) for x in vect], sep=the_sep)
    f.write(line)
    f.close()

def read_vect_to_float(f, the_sep = ','):
    r = csv.reader(f, delimter = the_sep, quotechar = '')
    line = r.next()
    vect = [float(x) for x in line]
    r.close()
    return vect

def read_mat_to_float(f, the_sep = ','):
    r = csv.reader(f, delimter = the_sep, quotechar = '')
    mat = []
    for line in r:
        vect = [float(x) for x in line]
        mat.append(vect)
    r.close()
    return mat

def read_vect_to_int(f, the_sep = ','):
    r = csv.reader(f, delimter = the_sep, quotechar = '')
    line = r.next()
    vect = [int(x) for x in line]
    r.close()
    return vect

def read_mat_to_int(f, the_sep = ','):
    r = csv.reader(f, delimter = the_sep, quotechar = '')
    mat = []
    for line in r:
        vect = [int(x) for x in line]
        mat.append(vect)
    r.close()
    return mat

def get_representative_atom(res):
    if 'CA' in res.child_dict.keys():
        return res['CA']
    elif 'CB' in res.child_dict.keys():
        return res['CB']
    else:
        print 'no CA or CB atom in residue'
            #pdb.set_trace()
        return res.child_list[0]
