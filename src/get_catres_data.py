import os

home = '../catres-fam/'

pdb_folders = os.listdir(home)
pdb_chains = []


pdb_files_f = open('../catres_pdbs', 'w')
master_f = open('../catres_sites','w')




for folder in pdb_folders:
    s = folder.split(':')
    pdb = s[0]
    chain = ''
    if len(s) == 2:
        chain = s[1]
    else:
        chain = '-1'
    pdb_files_f.write(pdb + '_' + chain + '\n')

    # open directory to get the residues
    info_file = home + folder + '/' + pdb + '.catres'
    f = open(info_file)
    for line in f:
        w = line.split('\t')[0].split(' ')
        temp_letter = w[0]
        pos = w[2]
        master_f.write(pdb + ',' + chain + ',' + str(pos) + '\n')


