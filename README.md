# ZIPcnv: accurate and efficient inference of copy number variations from shallow whole-genome sequencing

## Abstract
Shallow whole-genome sequencing (sWGS), a rapid and cost-effective sequencing technology, has gradually been widely adopted for CNV analyses. Howev er, with genome‑wide coverage of only 0.1–5X, sWGS data display a pronounced zero‑inflation phenomenon—a large fraction of loci has zero sequencing reads. Zeo‑inflation causes read counts to fluctuate by several‑fold between adjacent windows. As a result, random upward blips in coverage can be misinterpreted as copy‑number gains (false positives), and true deletions often become indistinguishable from pervasive zero‑coverage noise. In addition,  existing CNV detection tools developed for sWGS data often struggle to adapt across different CNV sizes. These combined effects severely constrain the accuracy of CNV inference. To address above challenges, we propose ZIPcnv, a novel CNV detection tool specifically designed for sWGS data. First , we apply a large sliding window to smooth the raw read depth signal, which transforms the original zero-inflated statistical characteristics into approximately normal distribution characteristics. We  then design a statistical process model that robustly detects persistent shifts under high background noise using a cumulative sum strategy, classifying genomic regions into candidate and non-candidate CNV regions. Finally, dynamic sliding windows are used for one-pass detection of CNVs of varying lengths, with window size adapting to the CNV region size. We evaluated the performance of ZIPcnv on simulated data and 190 real whole-genome sequencing samples. Experimental results show that ZIPcnv consistently outperforms currently popular CNV detection tools.


## Installation
Uncompress the installation package using the following commands:

```bash
cd /my/install/dir/
unzip /path/to/ZIPcnv.zip
```

**Requirements**

You'll need to install the following Python packages in order to run the code:

```bash
Python 3.11.9 
pysam 0.22.1
json 1.9.6
tqdm 4.67.1
pandas 2.2.3
numpy 2.2.3
scipy 1.15.2
```

Before starting the project, you need to configure the parameter file. For detailed instructions, refer to **my.config**. 


## Step 1: Data Preprocessing
The preprocessing consists of two steps: 1.data normalization and 2.baseline setting.

**Usage**
```bash
python3 data_processing.py [-config CONFIG]

commands:
-config [str]: Path to the configuration file.
```

Example:
```bash
python3 data_processing.py -config my.config
```

Note: Setting the baseline requires at least 50 normal samples; otherwise, a warning will be issued.

## Step 2: CNV detection using dynamic sliding windows 
ZIPcnv employs a CUSUM control chart-based model for CNV detection.

**Usage**
```bash
python3 zip-caller.py [-config CONFIG] [-o OUTPUT] [-w SLIDING_WINDOW_SIZE] [-k REFERENCE_VALUE]

commands:
-config [str]: Path to the configuration file.
-o [str]: Path to the output file.
-w [int]: Sliding window size (default: 3000).
-k [float]: Reference value for the allowed degree of deviation (default: 0.3).
```

Example:
```bash
python3 zip_caller.py -config my.config -o data/zipcall-output -w 3000 -k 0.3
```

ZIPcnv is an intermediate step of PGcnv. For a more detailed introduction to ZIPcnv, please refer to PGcnv (https://github.com/Nevermore233/PGcnv).
