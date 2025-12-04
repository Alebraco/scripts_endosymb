import pandas as pd
import os
import Bio.Phylo as Phylo
import time
def distance_matrix(group):
    '''Creates a matrix based on patristic distance
    args:
        group (str): endosymbiont_only, endosymbiont+relatives, relatives_only 
    output:
        dist_df (df): DataFrame with distances matrix
    '''
    
    dist_matrices = {}
    tree_dir = os.path.join(group,'dna_tree_results/')
    
    print(f'Starting distance matrix computation for {group}')
    print(f'Processing trees in {tree_dir}')

    species_list = os.listdir(tree_dir)
    species_number = len(species_list)
    print(f'Found {species_number} species directories')
    
    i = 1
    for sp in os.listdir(tree_dir):
        print(f'{i}/{species_number} Processing {sp}') 
        start_time = time.time()
        tree_path = os.path.join(tree_dir, sp, f'{sp}.treefile')
        i += 1
        if not os.path.exists(tree_path):
            print(f'No tree file found for {sp} at {tree_path}.')
            continue

        tree = Phylo.read(tree_path, "newick")
        terminals = tree.get_terminals()
        print(f'{len(terminals)} terminals were found.')

        dist_matrix = [[tree.distance(t1, t2) for t2 in terminals] for t1 in terminals]

        print(f'Done. Computation time was {time.time() - start_time}')
        ids = [t.name for t in terminals]
        dist_df = pd.DataFrame(dist_matrix, index = ids, columns = ids)
        dist_matrices[sp] = dist_df
    print(f'{len(dist_matrices)} matrices were computed.')
    return dist_matrices
