# antibodyConcat
The procedure are currently split into parts for easier modification.

## I_generateInputReads.py
As the name said, it selects the input reads from the original input files. 
You can define the parameters of your splicing by running this code.
This code will generate a folder named by parameters, which also contains the setting(parameters) file and the input reads.

## II_assembleFromReads.py
The start step of splicing, you can just feed the folder path that was generated by the "I" step.
The output is the initial splicing output waiting for further processing.
### 
