import pdb
import f as features
import new_new_objects as objects
import param
import wc
import global_stuff

# hardcode parameters for the experiment here for now.

#the_dict = {'pdb_name':'1asy', 'chain_letter':'A', 'edge_feature_list':[features.xW], 'node_feature_list':[features.vW, features.uW, features.wW], 'dist_cut_off':5}

#the_params = param.param(the_dict)
from parameters import the_params
pdb.set_trace()
ans = wc.get_stuff(objects.bmW, the_params, global_stuff.recalculate, False, False)

asdf = [ans[j][0] for j in range(len(ans))]
sum = 0
for i in asdf:
    sum += i*i
print sum / len(ans)
