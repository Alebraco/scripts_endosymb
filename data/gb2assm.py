#!/usr/bin/env python3

import os
from Bio import Entrez
import time
import json

dataset_path = 'files/endosymbiont_dict.json'
with open(dataset_path, 'r') as f:
    dataset_dict = json.load(f)
existing_gcfs = set(sum(dataset_dict.values(), []))

Entrez.email = 'asoneto@ncsu.edu'

related_sp = 'related_species/'
for sp in os.listdir(related_sp):
	print(f'Processing {sp}:')
	sp_path = os.path.join(related_sp, sp)
	for file in os.listdir(sp_path):
		print(f'\tProcessing file: {file}')
		if not file.endswith('_host.tsv') and not file.endswith('accns.txt'):
			file_path = os.path.join(sp_path, file)
			new_file_path = file_path.replace('.tsv', 'accns.txt')
			gb_to_pident = {}
			with open(file_path, 'r') as f:
				for line in f:
					if line.strip():
						gb_acc, pident = line.strip().split('\t')[:2]
						gb_to_pident[gb_acc] = pident

			gcf_to_pident = {}
			for gb, pid in gb_to_pident.items():
				try:
					handle = Entrez.elink(
						dbfrom='nuccore',
						db='assembly',
						id=gb
					)
					record = Entrez.read(handle)
					handle.close()
					time.sleep(0.3)

					if record and len(record) > 0 and 'LinkSetDb' in record[0]:
						for linkset in record[0]['LinkSetDb']:
							for link in linkset['Link']:
								gcf_id = link['Id']
								handle2 = Entrez.esummary(db='assembly', id=gcf_id)
								summary = Entrez.read(handle2)
								handle2.close()
								time.sleep(0.3)

								for doc in summary['DocumentSummarySet']['DocumentSummary']:
									assembly_acc = doc.get('AssemblyAccession')
									if assembly_acc and assembly_acc.startswith('GCF') and assembly_acc not in existing_gcfs:
										gcf_to_pident[assembly_acc] = pid
										print(f'\t >Processed {gb} with pident {pid}')
				except Exception as e:
					print(f'\t >Failed to process {gb}: {e}')
					continue
			with open(new_file_path, 'w') as f:
				for gcf, pid in gcf_to_pident.items():
					f.write(f'{gcf}\t{pid}\n')
