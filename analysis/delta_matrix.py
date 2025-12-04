#!/usr/bin/env python3
import os
import json
import pandas as pd
import numpy as np
from utils import files_dir

def delta_matrix(data, mat_type, save_to_file = False, filename = None):
    '''Create a delta matrix
    args:
        data(dict): nested dictionary {species: {id: mat_type}}
        save_to_file(boolean)
        filename(str)
    output:
        delta_matrices: delta matrix for each species specified in data
    '''
    delta_matrices = {}


    if save_to_file:
        if filename is None:
            filename = os.path.join(files_dir, f'dmatrix_{mat_type}.json')
        try:
            with open(filename, 'r') as save:
                delta_matrices_data = json.load(save)
            for sp, data in delta_matrices_data.items():
                delta_matrices[sp] = pd.DataFrame(data['matrix'], index=data['index'], columns=data['columns'])
            print('Delta matrices loaded from', filename)
            return delta_matrices

        except:
            print('No delta matrices available. Creating new matrices.')



    for sp, dataset in data.items():
        ids = list(dataset.keys())
        gc_values = [dataset[acc][mat_type] for acc in ids]

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
