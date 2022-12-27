# antibodyConcat
The procedure are currently split into parts for easier modification.

## Installation
Simply using git clone.

## Requirements
    Predfull
    Rapsearch
    msslash
We recommend to install Predfull and msslash by their github's guidance(git clone and install), and install the rapsearch by conda (e.g. conda install rapsearch)

## Sample usage
Assuming the Predfull and msslash are installed at the same path as this tool:

    python run.py -source avastin -score 0.8 -t 2 -kl 5 -ku 8 -more 0 -predfull_path PredFull/predfull.py -msslash_path msSLASH/bin/bruteforce
It may take few minutes, the final results will appeal on the generated folder named by parameters, which contains an HTML report and two concatenation result for heavy chain and light chain separately. 

