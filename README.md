PBPM
==============================

General protocol to generate patient by disrupted pathways matrices (PBPMs).

Project Organization
------------

    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── analysis           <- The results and files created by the code, cytoscape sessions, etc.
    │
    ├── data
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── environment.yml   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `conda env export --no-builds > env.yml`
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Logs, data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── src                <- Source code for use in this project.
    │   └── __init__.py    <- Makes src a Python module


### Description

This repo contains a Python module to develop a general protocol to be executed to generate patient by disrupted 
pathways matrices (also known as PBPMs).

### Interface

This protocol is based on a function-oriented program, whose functions are called by a main method and executed from
 the command line interface. Future devops plans include the possibility of building a web GUI.

This project contains a module called [PBPM.py](hhttps://github.com/mgmartinezl/Stalicla-PBPM/blob/master/src/PBPM.py), which contains
a set of functions integrated and called by the [main-PBPM.py](hhttps://github.com/mgmartinezl/Stalicla-PBPM/blob/master/src/main-PBPM.py) 
script. 

#### Positional parameters

* **inputFile:** it is mandatory to specify the absolute path to the input file containing the 
patient mutations.
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt```
* **pathwaysDirectory:**  it is mandatory to specify the absolute path to the directory that 
contains the pathway files with gene annotations.
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/```

#### Optional parameters

* **-pathway:** optional argument to extract only specific pathways. If no setting is provided,
all available pathways will be extracted by default.
    - Example 1: it will run only for pathway R-HSA-69620 \
    ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -pathway R-HSA-69620 ```
    - Example 2: it will run for all possible pathways \
    ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ ```    
    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one pathway, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of pathways separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-pathway R-HSA-69620,0051705 ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired pathways to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-pathway ~/PBPM/data/raw/filters/my-file-containing-pathways.txt``` \

* **-gene:** optional argument to extract only specific genes. If no setting is provided,
all available genes will be extracted by default.
    - Example: it will run only for gene CTR9 \
    ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -gene CTR9 ```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one gene, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of genes separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-gene CTR9,NOCL2 ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired genes to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-gene ~/PBPM/data/raw/filters/my-file-containing-genes.txt``` 

* **-patient:** optional argument to extract only specific patient IDs. If no setting is provided,
all available patients will be extracted by default.
    - Example: it will run only for patient with ID 1 \
    ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -patient Patient_1 ```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one patient ID, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of patients separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-patient Patient_X,Patient_Y ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired patient IDs to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-patient ~/PBPM/data/raw/filters/my-file-containing-patients.txt``` 

* **-mutation:** optional argument to extract only specific mutations. If no setting is provided,
all available consequences will be processed by default.
    - Example: it will run only for mutations of type 'missense_variant' \
    ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -mutation missense_variant```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one mutation, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of mutations separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-mutation missense_variant,Intron ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired mutations to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-mutation ~/PBPM/data/raw/filters/my-file-containing-mutations.txt``` 

* **-pli_gt:** optional argument to filter records with values greater than or equal to a specified pLI threshold.
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -pli_gt 0.6 ```
   
* **-af_lt:** optional argument to filter records with values less than a max_control_AF threshold.
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -af_lt 0.3 ```

* **-pph2:** optional argument to filter records for qualifiers of pph2 predictions. Available options for this
parameter are: *benign*, *possibly damaging*, and *probably damaging*. Note that qualifiers must be enclosed with 
quotes as they hold blank spaces.
    - Example:  ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -pph2 "possibly damaging" ```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; To filter more than one pph2 qualifier, you must separate their labels with comma and without spaces while keeping the quotes: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ``` -pph2 "probably damaging,possibly damaging" ```

* **-mpc_gt:** optional argument to filter records with values greater than or equal to a specified MPC threshold.
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -mpc_gt 0.7 ```

* **-adj:** optional argument to filter records by specific value(s) of adjusted consequence. Available options for this
parameter include: *PTV*, *Missense3*, *Missense*, etc. 
    - Example:  ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -adj PTV ```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; To filter more than one pph2 qualifier, you must separate their labels with comma: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ``` -adj Missense3,Missense ```

* **-matrix:** optional argument to select the type of matrix to be extracted. It is only possible to generate one matrix at a time.
Available options are: *binary(b)*, *numerical(n)*, *normalized(nn)*. If not set, the default value of this argument will return a binary matrix.
    - Example to explicitly generate a binary matrix : \
     ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -matrix b ```
    - Example to generate a numerical matrix (not normalized): \
      ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -matrix n ```
    - Example to generate a normalized numerical matrix: \
     ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -matrix nn ```

* **-r:** set this parameter to yes (y) to download a copy of the base matrix (raw data) generated to build the PBP matrix. 
By default, this parameter is set to no (n).
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways/ -r y ```

* **-append:** when this parameter is specified, the base matrix of the current session will be appended to the file provided. 
Thus, it will not be necessary to explicitly set the -r parameter to 'yes', unless a separate copy of the base matrix is desired.
    - Example: ```$ python3 main-PBPM.py ~/PBPM/data/raw/original-mutations-file.txt ~/PBPM/data/raw/pathways// -append path-to-my-historic-file.txt ```
  

**Note:** help() is available for all the parameters via the command line. 

### Input files

A sample of the file required by the first positional argument (the input file containing mutations and patients)
can be found in the folder **data/raw**. For direct access, click [here.]() 

Similarly, in the **data/raw/pathways** folder, a sample of a directory containing pathways can be found, which is a mandatory
parameter for _pathwaysDirectory_ entry. See it directly [here.]()

Additional example text files to filter can be found in the folder **data/raw/filters**
[genes](), 
[patients](), 
[pathways]() and
[mutations]() 
can be found in the folder as well.

### Running tests

In the folders [references]() and 
[data/processed](), three different running
examples can be found. Each of them generates a log containing the parameters set to run the
program, as well as the desired output.

## How to run this script

The scripts PBPM.py, and main-PBPM.py are written in Python 3, 
which uses up-to-date libraries for this version as well.
 
To run the main module in a linux environment, simply call the script and the
arguments it needs:

```python3 main-PBPM.py inputFile pathwaysDirectory -pathway -gene -patient -mutation -pli_gt -af_lt -pph2 -mpc_gt adj -matrix -r -append```

For any additional information, contact me: 

*Gabriela Martinez* <br>
*airamgabriela17@gmail.com*

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
