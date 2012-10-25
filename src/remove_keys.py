import constants

import subprocess

wrappers = ['bwW', 'caW', 'cbW', 'ccW', 'cdW']

for wrapper in wrappers:
    subprocess.call('rm' + ' ' + constants.BIN_FOLDER + '*' + wrapper + '*', shell=True)

subprocess.call('rm ../new_results/2_0', shell=True)
