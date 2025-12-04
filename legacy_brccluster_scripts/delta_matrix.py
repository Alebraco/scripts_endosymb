#!/usr/bin/env python3
import os
import json
import pandas as pd
import numpy as np

def delta_matrix(data, type, save_to_file = False, filename = f'dmatrix_{type}.json'):
    '''Create a delta matrix
    args:
        data(dict): nested dictionary {species: {id: type}}
        save_to_file(boolean)
        filename(str)
    output:
        delta_matrices: delta matrix for each species specified in data
    '''
    delta_matrices = {}

    if save_to_file and os.path.isfile(filename):
        try:
            with open(filename, 'r') as save:
                delta_matrices_data = json.load(save)
            for sp, data in delta_matrices_data.items():
                delta_matrices[sp] = pd.DataFrame(data['matrix'], index=data['index'], columns=data['columns'])
            print('Delta matrices loaded from', filename)
            return delta_matrices

        except:
            print('No delta matrices available. Creating new matrices.')



    for sp, data in data.items():
        ids = list(data.keys())
        gc_values = [data[id][type] for id in ids]

        dims = len(ids)
        if dims < 2: 
            print(f'Only one genome is available for {sp}')
            continue
        delta_matrix = np.zeros((dims,dims))

        for i in range(dims):
            for j in range(dims):
                delta_matrix[i][j] = abs(gc_values[i] - gc_values[j])

        delta_df = pd.DataFrame(delta_matrix, index = ids, columns = ids)
        delta_matrices[sp] = delta_df

    if save_to_file:
        delta_matrices_data = {}
        for sp, df in delta_matrices.items():
            delta_matrices_data[sp] = {
                'matrix': df.values.tolist(),
                'index': df.index.tolist(),
                'columns': df.columns.tolist()
            }
        with open(filename, 'w') as save:
            json.dump(delta_matrices_data, save)
        print('Delta matrices saved to', filename)

    return delta_matrices
