#!/usr/bin/env python3

import os
import pandas as pd
from utils import files_dir

def load_IS_sequences(file):
    IS_data = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                IS_data.append(line.strip())


if os.isfile(os.path.join(files_dir, 'IS_sequences.faa'):
    
