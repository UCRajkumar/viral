# HPV Detection

## Authors: Utkrisht Rajkumar, Dongxia Wu

### Description
In this project, We use a variational autoencoder (VAE) to detect HPV reads in human tumour genome data. We train our VAE on HPV sequences and test on human tumour genome data. We use a training set with only HPV genomes because real world datasets contain reads from rapidly evolving viral, bacterial, fungal, and human genomes. It is nearly impossible to create a fully inclusive training set that is representative of any future tumor genome dataset. Therefore, we aim to create a model that performs openset recognition. 

## Table of contents

- [Requirements](#requirements)
- [Preprorcessing](#preprocessing)
- [How to run](#run)
- [Description of files](#description)

## Requirements <a name="requirements"></a>
* tensorflow 2.3.1
* numpy 1.18.5
* pandas 1.1.3
* keras 2.4.0
* sklearn 0.23.2
* Bio 1.78

## Preprocessing <a name="preprocessing"></a>
Please contact urajkuma [at] eng.ucsd.edu for data.

## How to run  <a name="run"></a>

Run the following scripts to do openset recognition test:
```
temporal_vanilla_vae.ipynb, beta_vae.ipynb, sequential_vanilla_vae.ipynb, Evaluation.ipynb
```

## Description of files  <a name="description"></a>
 

file name | Description 
--- | ---
Evaluation.ipynb | Script to plot ROC and calcualate AUC as a metric for comparison.
beta_vae_ELBOW.ipynb | Script for beta VAE model.
temporal_vanilla_vae.ipynb | Script for vanilla VAE model.
sequential_vanilla_vae.ipynb | Script for sequential VAE model.
data_preprocess.ipynb | Script to visualize reads from viral, bacterial, fungal, and human genomes.
utils.py | Contains read function, label assignment function, and one-hot encoding function.




