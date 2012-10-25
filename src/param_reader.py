from param import paramm
import constants
import global_stuff


def read_param(file_location):

    # read folder_name
    f = open(file_location)
    the_params = param.param()
    hp_values = param.param()
    folder_name = f.readline().strip()

    

    for line in f:
        s = line.strip().split(',')
        if s[0] != 'hp':
            the_type = s[1]
            the_name = s[0]
            if the_type == 'f':
                the_params.set_param(the_name, float(s[2]))
            elif the_type == 'i':
                the_params.set_param(the_name, int(s[2]))

    the_params.set_param('hp', hp_values)
    return folder_name, the_params

def read_info(file_location):

    # read folder_name
    f = open(file_location)
    the_params = param.param()
    hp_values = param.param()
    folder_name = f.readline().strip()

    

    for line in f:
        s = line.strip().split(',')
        if s[0] == 'hp':
            to_add = []
            the_type = s[2]
            the_name = s[1]
            for i in range(2,len(s)):
                if the_type == 'f':
                    to_add.append(float(s[i]))
                elif the_type == 'i':
                    to_add.append(int(s[i]))
            hp_values.set_param(the_name, to_add)
        else:
            the_type = s[1]
            the_name = s[0]
            if the_type == 'f':
                the_params.set_param(the_name, float(s[2]))
            elif the_type == 'i':
                the_params.set_param(the_name, int(s[2]))

    the_params.set_param('hp', hp_values)
