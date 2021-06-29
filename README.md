# DeconToVCF

## Introduction

This is a tool to convert multiple DeCon .txt file outputs generated in pipelines into a single genomic VCF format file that can then be used for annotation software and further analysis tools.

## Download the directory:
```
git clone git@github.com:AWGL/DeconToVCF.git
```

## Requirements

The python script needs to be copied to and run from within pipeline post-processing folder.

The post-processing folder should contain a "Results" folder which contains:
* "ped" folder containing "<runid>.ped" file
* "cnv_svs/raw_data" folder containing multiple .txt reports of genomic variant data outputted from DeCon


## To run the programme

```
python DeconToVCF
```

## Additional info

* Output will be a file within the directory named "<runid>_full_CNV.cnv.vcf"
