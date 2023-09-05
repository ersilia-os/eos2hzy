import requests
from standardiser import standardise
from rdkit import Chem
import urllib

def run_chemed(origin_smiles: str, similarity: float = 0.7):
    """Function adapted from Andrew White's Exmol"""

    url_smiles = urllib.parse.quote(origin_smiles)

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{url_smiles}/property/CanonicalSMILES/JSON"
    try:
        reply = requests.get(
            url,
            params={"Threshold": int(similarity *100), "MaxRecords": 100},
            headers={"accept": "text/json"},
            timeout=10,
        )
    except requests.exceptions.Timeout:
        print("Pubchem seems to be down right now")
        return []
    try:
        data = reply.json()
        if "PropertyTable" not in data or "Properties" not in data["PropertyTable"]:
            return []
        smiles = [d["CanonicalSMILES"] for d in data["PropertyTable"]["Properties"]]
        smiles = list(set(smiles))
        return smiles
    except:
        return []
    


class PubChemSampler(object):
    def __init__(self):
        pass

    def _sample(self, smiles):
        smiles = run_chemed(origin_smiles=smiles)
        std_smiles = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                mol = standardise.run(mol)
            except:
                continue
            std_smiles += [Chem.MolToSmiles(mol)]
        return std_smiles

    def sample(self, smiles):
        sampled_smiles = []
        sampled_smiles += self._sample(smiles)
        sampled_smiles = list(set(sampled_smiles))
        return sampled_smiles