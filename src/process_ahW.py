import glob

import global_stuff

results_folder = global_stuff.RESULTS_FOLDER
files = glob.glob(results_folder + 'mfdW-ahW*')

for a in files:
    f = open(a,'r')
    s = f.readlines()[-1]
    print s.strip(), a
    
