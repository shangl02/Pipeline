# scEpiLock
Pytorch implementation for **scEpiLock: a weakly supervised learning scheme for cis-regulatory element localization and variant impact quantification for single-cell epigenetic data**

## Setup
Clone the repository.
```
git clone https://github.com/yanwengong/scEpiLock/
```
Set up the conda envrionment. Note the spec-file.txt is for Linux platform. 
```
conda create --name scEpilock --file scEpilock_spec-file.txt
```
Activate conda environment
```
conda activate scEpilock
```

## Data
Please download the data [here](https://drive.google.com/file/d/1hS1bBdlAxZjYOAqgenUe62oX79qIDCYD/view?usp=sharing). 

## Multi-label Deep Learning Module 
This is to train and test the multi-label deep learning module, which predicts chromatin accessibility for individual cell type.
### Configuration 
Configuration json files are stored at scEpiLock/deep_learning_module/config/. For instance, you can modify the brain_config.json based on the input location and desired result location. The following feilds need to be modified.
- "pos_forward_path": "/absolute/path/to/brain/brain_atac_hg38_final.fa",
  - this is the location for the input fasta of the scATAC-seq peak 
- "neg_path": "/absolute/path/to/brain/neg_noENCODE_noBrain.fa",
  - this is the location for the negative control pasta file
- "label_path": "/absolute/path/to/brain/label_one_hot.txt"
  - this is the location for the scATAC-seq peak label file
- "model_path": "/absolute/path/to/brain/deep_learning_res/model.pt",
  - where the trained model should be saved
- "output_evaluation_data_path": "/absolute/path/to/brain/deep_learning_res/"
  - where the result should be saved 

### Train
```
python scEpiLock/deep_learning_module/main.py  -e 'train' -s 'main' -c 'scEpiLock/deep_learning_module/config/brain_config.json'
```

### Test
```
python scEpiLock/deep_learning_module/main.py  -e 'test' -s 'main' -c 'scEpiLock/deep_learning_module/config/PBMC_config.json'
```

## Object Detection Module 
This is a weakly supervised, Grad-CAM module that precisely locates the cis-regulatory element boundires. 

### Configuration 
Configuration json files are stored at scEpiLock/grad_cam_module/config/. For instance, you can modify the brain_config.json based on the input location and desired result location. The following feilds need to be modified.
- "pos_forward_path": "/absolute/path/to/brain/brain_atac_hg38_final.fa",
  - this is the location for the input fasta of the scATAC-seq peak 
- "model_path": "/absolute/path/to/brain/deep_learning_res/model.pt",
  - where the previously trained model is saved
- "result_path": "/absolute/path/to/brain/cam_res"
  - where the object-detection module result should be saved

### Run Object Detection Module
```
python scEpiLock/object_detection_module/main.py config/brain_config.json
```

## Variant Impact Quantification Module
This module quantifies the non-coding variant impact on chromatin accessibility for each cell type.

### Configuration 
Configuration json files are stored at scEpiLock/variant_impact/config/. For instance, you can modify the brain_config.json based on the input location and desired result location. The following feilds need to be modified.
- "peak_snp_bed": "/absolute/path/to/brain/intersect_peak_snp_loc_all_sort.bed",
- "ref_fasta": "/absolute/path/to/brain/hg38.fa",
- "model_weight_path": "/absolute/path/to/brain/deep_learning_res/model.pt",
  - where the previously trained model is saved
- "result_path": "/absolute/path/to/brain/variant_impact/all_snp_20210819",
  - where the variant impact result should be saved
  
### Run Variant Impact Quantification Module
```
python scEpiLock/variant_impact/main.py scEpiLock/variant_impact/config/brain_config.json
```
