# PubChem Molecular Sampler

A simple sampler of the PubChem database using their API. It looks for similar molecules to the input molecule and returns a list of 100 molecules by default. This model has been developed by Ersilia and posts queries to an online server.

## Identifiers

* EOS model ID: `eos2hzy`
* Slug: `pubchem-sampler`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Similarity`
* Output: `Compound`
* Output Type: `String`
* Output Shape: `List`
* Interpretation: 100 nearest molecules in PubChem

## References

* [Publication](https://academic.oup.com/nar/article/51/D1/D1373/6777787)
* [Source Code](https://github.com/ersilia-os/chem-sampler/blob/main/chemsampler/samplers/pubchem/sampler.py)
* Ersilia contributor: [GemmaTuron](https://github.com/GemmaTuron)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos2hzy)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos2hzy.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos2hzy) (AMD64)

## Citation

If you use this model, please cite the [original authors](https://academic.oup.com/nar/article/51/D1/D1373/6777787) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a GPL-3.0 license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!