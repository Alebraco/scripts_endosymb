#!/usr/bin/env python3
import os
import json
maindir = 'endosymb_only/'
checks = ['endosymb+relatives/', 'relatives_only/']
file = 'exclude.json'

with open(file, 'r') as f:
    exclude = json.load(f)

ids = [id for slist in exclude.values() for id in slist]

for i in ids:
    for root,dirs,files in os.walk(maindir):
        for file in files:
            if file.startswith(i):
                path = os.path.join(root,file)
                os.remove(path)