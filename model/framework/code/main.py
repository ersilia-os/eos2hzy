# imports
import os
import csv
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from chemsampler.sampler import ChemSampler
# import argparse
import pandas as pd
import csv

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]
    
# run model
sampler = ChemSampler()

# in sampled_smiles_dict, the keys will be the input smiles and the values will be a list of similar smiles
sampled_smiles_dict = sampler.sample(
        smiles_list=smiles_list,
        Sampler='PubChemSampler',
        num_samples=100,
        sim_ub=0.9, # these sim_ub and sim_lb values are the recommended values from the chem-sampler README
        sim_lb=0.3,
        distribution='ramp',
    )

if sampler == "all":
    # result is a dict where keys are name of the sampler and values are the sampled smiles
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(sampled_smiles_dict)):
            key = list(sampled_smiles_dict.keys())[i]
            # writer.writerow([key]) # it will write a single list, if want to retrieve by sampler, uncomment this
            for j in range(len(sampled_smiles_dict[key])):
                writer.writerow([sampled_smiles_dict[key][j]])
else:
    max_similarities = max(len(sim_list) for sim_list in sampled_smiles_dict.values())  # find the number of column labels we need for the .csv file
    header = ["sim-{0}".format(i+1) for i in range(max_similarities)]
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for key in sampled_smiles_dict: 
            writer.writerow(sampled_smiles_dict[key])   # sampled_smiles_dict[key] is a list of the similar molecules 

#check input and output have the same length
# input_len = len(smiles_list)
# output_len = len(outputs)
# assert input_len == output_len
