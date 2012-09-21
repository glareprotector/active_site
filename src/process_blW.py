import glob

results_folder = '../results/'
files = glob.glob(results_folder + 'blW*')

for a in files:
    f = open(a,'r')
    s = f.readlines()[0]
    print s.strip(), a
    
