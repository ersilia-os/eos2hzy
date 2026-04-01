import time
import requests
from standardiser import standardise
from rdkit import Chem
import urllib
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

def calculate_similarity(ref_mol, mol_list):
    ref_fp = AllChem.GetMorganFingerprint(ref_mol, 2)
    mol_fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in mol_list]
    similarities = [
        DataStructs.TanimotoSimilarity(ref_fp, mol_fp) for mol_fp in mol_fps
    ]
    return similarities


def sort_molecules_by_similarity(ref_mol, mol_list):
    similarities = calculate_similarity(ref_mol, mol_list)
    paired = list(zip(mol_list, similarities))
    paired.sort(key=lambda x: (-x[1], Chem.MolToSmiles(x[0], canonical=True)))
    sorted_mols = [mol for mol, sim in paired]
    return sorted_mols


class PubChemSampler(object):
    def __init__(self):
        pass

    def _run_chemed(self, origin_smiles: str, similarity: float = 0.7):
        """Function adapted from Andrew White's Exmol"""
        url_smiles = urllib.parse.quote(origin_smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{url_smiles}/property/CanonicalSMILES/JSON"
        for attempt in range(10):
            try:
                reply = requests.get(
                    url,
                    params={"Threshold": int(similarity*100), "MaxRecords": 1000},
                    headers={"accept": "text/json"},
                    timeout=10,
                )
            except requests.exceptions.RequestException as e:
                # Network-level failure (timeout, DNS error, connection refused, etc.)
                print(f"PubChem not accessible for compound '{origin_smiles}' (attempt {attempt+1}/10): {e}. Waiting 30s...")
                time.sleep(30)
                continue

            if reply.status_code == 429:
                # Rate limit exceeded — PubChem allows ~5 requests/second; back off and retry
                print(f"PubChem rate limit hit for compound '{origin_smiles}' (attempt {attempt+1}/10). Waiting 60s...")
                time.sleep(60)
                continue

            if reply.status_code != 200:
                # Unexpected HTTP error (e.g. 400 bad request, 500 server error) — no point retrying
                print(f"PubChem returned HTTP {reply.status_code} for compound '{origin_smiles}'. Skipping.")
                return []

            # Successful response — parse and return
            try:
                data = reply.json()
                if "PropertyTable" not in data or "Properties" not in data["PropertyTable"]:
                    # Valid response but no similar compounds found
                    return []
                smiles = [d["ConnectivitySMILES"] for d in data["PropertyTable"]["Properties"]]
                smiles = list(set(smiles))
                return smiles
            except Exception as e:
                # Response was not valid JSON or had unexpected structure
                print(f"Failed to parse PubChem response for compound '{origin_smiles}': {e}. Skipping.")
                return []

        # All 10 attempts exhausted without a successful response
        print(f"PubChem unreachable after 10 attempts for compound '{origin_smiles}'. Giving up.")
        return []
    

    def _sample(self, origin_smiles):
        sampled_smiles = self._run_chemed(origin_smiles=origin_smiles)
        std_smiles = []
        for smi in sampled_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                mol = standardise.run(mol)
            except:
                continue
            std_smiles += [Chem.MolToSmiles(mol)]
        return std_smiles

    def _sort_smiles(self, origin_smiles, sampled_smiles):
        ref_mol = Chem.MolFromSmiles(origin_smiles)
        if len(sampled_smiles) == 0:
            return []
        mol_list = []
        for smi in sampled_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol_list += [mol]
        sorted_mols = sort_molecules_by_similarity(ref_mol, mol_list)
        sorted_smiles = []
        for m in sorted_mols:
            if m is None:
                continue
            sorted_smiles += [Chem.MolToSmiles(m)]
        return sorted_smiles


    def sample(self, origin_smiles):
        sampled_smiles = []
        sampled_smiles += self._sample(origin_smiles)
        sampled_smiles = list(set(sampled_smiles))
        sampled_smiles_ordered = self._sort_smiles(origin_smiles, sampled_smiles)
        sampled_smiles_selected = sampled_smiles_ordered[:100]
        if len(sampled_smiles_selected) < 100:
            sampled_smiles_selected.extend([None] * (100 - len(sampled_smiles_selected)))
        return sampled_smiles_selected