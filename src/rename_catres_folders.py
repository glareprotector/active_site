import os
import objects
import subprocess
import global_stuff

home = '../catres-fam/'

pdb_folders = os.listdir(home)

for folder in pdb_folders:
    s = folder.split(':')
    pdb = s[0]
    print s
    if s[1] == '':
        try:
            params = objects.param({'pdb_name':s[0]})
            structure = global_stuff.the_obj_manager.get_variable(objects.pdb_obj_wrapper(params))
            letter = structure[0].child_dict.keys()[0]
            print letter
            #rename
            subprocess.call(['mv', home + folder, home + folder + str(letter)])
        except:
            print 'fail'
