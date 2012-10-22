import os
import wc
import param
import new_new_objects as objects
import subprocess
import Bio.PDB
import pdb

home = '../data/catres-fischer/'

pdb_list_file = 'cf.pl'

pdb_folders = os.listdir(home)

f = open(pdb_list_file, 'w')

new_folders = []



for folder in pdb_folders:

    s = folder.strip().split('_')


    pdb_name = s[0]
    params = param.param({'p':pdb_name})

    while 1:
        try:
            g = wc.get_stuff(objects.fW, params, False, False, False)
        except Exception, err:
            print err
            import time
            time.sleep(20)
        else:
            break

    structure = Bio.PDB.PDBParser().get_structure(params.get_param('p'), g)
    if s[1] == '':
        letter = structure[0].child_dict.keys()[0]
    else:
        letter = s[1]
    params.set_param('c', letter)
    if s[2] == '':
        
        min_pos, max_pos = wc.get_stuff(objects.cgW, params, True, False, False)
        params.set_param('st', min_pos)
        params.set_param('en', max_pos)

        aa_to_pos = wc.get_stuff(objects.aW, params, False, True, False)
        min_pos = aa_to_pos[0]
        max_pos = aa_to_pos[-1]
    else:

        #min_pos = int(s[2].split('-')[0])
        #max_pos = int(s[2].split('-')[1])
        min_pos = int(s[2])
        max_pos = int(s[3])
        params.set_param('st', min_pos)
        params.set_param('en', max_pos)
        aa_to_pos = wc.get_stuff(objects.aW, params, False, True, False)

    pdb_list_line = pdb_name + ',' + letter + ',' + str(min_pos) + ',' + str(max_pos)
    f.write(pdb_list_line + '\n')
    new_folders.append(pdb_name + '_' + letter + '_' + str(min_pos) + '_' + str(max_pos))
    
# also rename the folder
"""
for i in range(len(new_folders)):
    folder = pdb_folders[i]
    new_folder = new_folders[i]
    subprocess.call(['mv', home+folder, home+new_folder])
    
"""
