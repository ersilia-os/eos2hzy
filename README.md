# PubChem Molecular Sampler

Returns 100 compounds from PubChem similar to a query, obtained through PubChem's two-dimensional fast similarity service. PubChem aggregates deposited chemical records at a scale far beyond curated bioactivity databases, so neighbours may include screening compounds, natural products and reagents rather than only drug-like chemistry. The request is issued to the public API at run time, meaning results follow the live database and depend on network availability rather than a bundled snapshot.

This model was incorporated on 2023-08-10.Last packaged on 2026-04-27.

## Information
### Identifiers
- **Ersilia Identifier:** `eos2hzy`
- **Slug:** `pubchem-sampler`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Similarity search`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Similarity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `100`
- **Output Consistency:** `Fixed`
- **Interpretation:** List of 100 compounds from PubChem structurally similar to the query molecule.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smiles_00 | string |  | Compound index 0 queried with the PubChem API |
| smiles_01 | string |  | Compound index 1 queried with the PubChem API |
| smiles_02 | string |  | Compound index 2 queried with the PubChem API |
| smiles_03 | string |  | Compound index 3 queried with the PubChem API |
| smiles_04 | string |  | Compound index 4 queried with the PubChem API |
| smiles_05 | string |  | Compound index 5 queried with the PubChem API |
| smiles_06 | string |  | Compound index 6 queried with the PubChem API |
| smiles_07 | string |  | Compound index 7 queried with the PubChem API |
| smiles_08 | string |  | Compound index 8 queried with the PubChem API |
| smiles_09 | string |  | Compound index 9 queried with the PubChem API |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Online`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos2hzy](https://hub.docker.com/r/ersiliaos/eos2hzy)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos2hzy.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos2hzy.zip)

### Resource Consumption
- **Model Size (Mb):** `1`
- **Environment Size (Mb):** `444`
- **Image Size (Mb):** `494.25`

**Computational Performance (seconds):**
- 10 inputs: `41.41`
- 100 inputs: `489.32`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/ersilia-os/chem-sampler](https://github.com/ersilia-os/chem-sampler)
- **Publication**: [https://doi.org/10.1093/nar/gkac956](https://doi.org/10.1093/nar/gkac956)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2023`
- **Ersilia Contributor:** [GemmaTuron](https://github.com/GemmaTuron)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [GPL-3.0-or-later](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos2hzy
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos2hzy
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
