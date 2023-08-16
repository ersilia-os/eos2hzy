# imports
import os
import csv
import sys
from chemsampler.sampler import ChemSampler
# import argparse
import pandas as pd
import csv


# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# my model
# def my_model(smiles_list):
#     return [MolWt(Chem.MolFromSmiles(smi)) for smi in smiles_list]


# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model
sampler = ChemSampler()
sampled_smiles = sampler.sample(
        smiles_list=smiles_list,
        Sampler='PubChemSampler',
        num_samples=100,
        sim_ub=0.9,
        sim_lb=0.3,
        distribution='ramp',
    )

if sampler == "all":
    # result is a dict where keys are name of the sampler and values are the sampled smiles
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(sampled_smiles)):
            key = list(sampled_smiles.keys())[i]
            # writer.writerow([key]) # it will write a single list, if want to retrieve by sampler, uncomment this
            for j in range(len(sampled_smiles[key])):
                writer.writerow([sampled_smiles[key][j]])
else:
    # result is a list of sampled smiles 
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(sampled_smiles)):
            writer.writerow([sampled_smiles[i]])
# outputs = my_model(smiles_list)

#check input and output have the same lenght
# input_len = len(smiles_list)
# output_len = len(outputs)
# assert input_len == output_len

# write output in a .csv file
# with open(output_file, "w") as f:
#     writer = csv.writer(f)
#     writer.writerow(["value"])  # header
#     for o in outputs:
#         writer.writerow([o])
