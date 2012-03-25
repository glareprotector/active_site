
#from objects import param
#from objects import pdb_chain_wrapper

from global_stuff import *
from objects import *
import subprocess
import os

print 'asdf'
print os.environ['PATH']
subprocess.Popen('echo $SHELL',shell=True)
subprocess.Popen('echo $PATH',shell=True, executable='/bin/bash')

inherited=param({'pdb_name': '1asy', 'chain_letter': 'A', 'evalue':1e-10})
c = pdb_chain_blast_results_file_wrapper(inherited, True)
#c=pdb_chain_seq_file_wrapper(inherited)
#item=the_obj_manager.get_variable(c)
item = the_file_manager.get_file_handle(c)
print item
#the_file_manager.get_file_handle(c)
print 'asdf'
