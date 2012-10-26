import os
import os.path


old_string = 'node_feature_list'
new_string = 'n'


folders = [x for x in os.listdir('.') if os.path.isdir(x)]
import pdb
import shutil
import re
for folder in folders:

    folder_files = [folder + '/' + x for x in os.listdir(folder)]
    for f in folder_files:
        new = f
        new = re.sub(r'node_feature_list', 'n', new)
        new = re.sub(r'edge_feature_list', 'e', new)
        new = re.sub(r'dist_cut_off', 'co', new)
        if f != new:
            shutil.copy2(f, new)
