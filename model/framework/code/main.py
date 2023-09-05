# imports
import os
import csv
import sys
from sampler import PubChemSampler


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
sampler = PubChemSampler()
outputs = []
for smi in smiles_list:
    o = sampler.sample(smi)
    outputs += [o]

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow([f"smi_{i}" for i in range(100)])  # header
    for o in outputs:
        writer.writerow(o)