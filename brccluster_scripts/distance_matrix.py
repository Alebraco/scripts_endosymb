import pandas as pd
import os
import Bio.Phylo as Phylo

def distance_matrix(group):
    '''Creates a matrix based on patristic distance
    args:
        group (str): endosymbiont_only, endosymbiont+relatives, relatives_only 
    output:
        dist_df (df): DataFrame with distances matrix
    '''
    
    dist_matrices = {}
    tree_dir = os.path.join(group,'dna_tree_results/')
    for sp in os.listdir(tree_dir):
        tree_path = os.path.join(tree_dir, sp, f'{sp}.treefile')

        if not os.path.exists(tree_path):
            print(f'No tree file found for {sp} at {tree_path}.')
            continue

        tree = Phylo.read(tree_path, "newick")
        terminals = tree.get_terminals()
        dist_matrix = [[tree.distance(t1, t2) for t2 in terminals] for t1 in terminals]
        ids = [t.name for t in terminals]

        dist_df = pd.DataFrame(dist_matrix, index = ids, columns = ids)
        dist_matrices[sp] = dist_df
        
    return dist_matrices