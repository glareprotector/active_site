from global_stuff import get_object
from objects import pdb_chain_pairwise_distance_obj_wrapper, pdb_chain_site_sorted_distances_obj_wrapper

# for now assume can retrieve scores and true states using wrapper
def get_performance(pdb_name, chain_letter, cutoff):
    param_dict = {'pdb_name':pdb_name, 'chain_letter':chain_letter}
    params = param(param_dict)
    sorted_distances, sorted_indicies = get_object(pdb_chain_site_sorted_distances_obj_wrapper, params)
    scores = get_object(pdb_chain_site_scores, params)
    true_states = get_object(pdb_chain_true_states, params)
    
    num_sites = len(sorted_distances)
    num_active_sites = sum(true_states)

    # determine how many closest sites to look at before quitting
    num_look = active_sites
    
    # find the indicies of the top ranked sites by sorting scores
    temp = [[i,scores(i)] for i in range(num_sites)]
    sorted_temp = sorted(temp, key = lambda x: x[1])
    best_indicies = [sorted_temp[i,0] for i in range(num_look)]


    num_found = 0

    for i in range(num_sites):
        if true_states[i] == 1:
            close_indicies = [sorted_indicies[i] for i in range(num_sites) if sorted_distances[i] < cutoff]
            if len(set(close_indicies) & set(best_indicies)) > 0:
                num_found = num_found + 1

    return num_found, num_active
