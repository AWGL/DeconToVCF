# DeconToVCF

## Introduction

This is a tool to convert multiple DeCon .txt file outputs generated in pipelines into a single genomic VCF format file that can then be used for annotation software and further analysis tools.

## Download the directory:
```
git clone git@github.com:AWGL/DeconToVCF.git
```

## Requirements

The python script needs to be run with arguments

* -d = path to decon output directory. The directory should be the folder containing all of the decon outputs "raw_data".
* -p = path to ped file.
* -o = path and filename of choice for the output. The programme will not create directories that do not already exist.


## To run the programme

* create environment from .yaml file and activate
```
conda env create --file DeconToVCF.yaml
conda activate DeconToVCF
```

* run programme with arguments
```
python DeconToVCF -d <path-to-decon-output-directory> -p <path-to-ped-file> -o <path/filename-for-output-file>
```
  eg:
```
python DeconToVCF -d /data/runid/post-processing/results/cnv_svs/raw_data/ -p post-processing/results/ped/file.ped -o post-processing/results/DeconToVCF_output.vcf
```

