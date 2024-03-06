# Bioinfo Project 2023

This repository contains materials for the bioinfomatics project (2023 L3 INFO Paris-Saclay).
The goal of this project is to identify biomarkers associated with ALS disease. 
To do this, you have access to RNA-Seq sequencing data from post-mortem brain cortex biopsies from individuals who had ALS and others who did not.

The data for this project comes from the study "Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon Activation, Oxidative Stress, and Activated Glia" by Tam et al. The full study can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6866666/).

## Introduction

This README will guide you throughout the project and your analyses. All steps below must be completed. Unless explicitly stated otherwise in the document (E.g., the mandaotry instructions), you have full freedom regarding the different options available to accomplish all tasks, especially concerning the coding.
This document will be updated at each session. 
Refer to the update date above.

You will work in pairs, from the begining. 
You are allowed to discuss with other pairs, but no direct code sharing (copy/paste etc). 
During all the project, note precisely your personnal contribution and the contribution of your pair.
You be asked for this during the last session (~15 min oral presentation, see below). 

During each of the first sessions, you will be given lessons/introductions on notions you may not have yet. 
Once this lessons start, stop all of your activities and follow the lessons, even if you are familar with the concept.

At the end of your project, you must send by email (phd.rinaudo@gmail.com):
- All your code,
- A report containing the result of your analyses,
- A text file (.csv, .txt...) contaning the ordered list of your top 100 genes (first = better).

The email must have the subject: "Projet_bioinfo_2023_nom1_prenom1_nom2_prenom2" and contains an archive named: "nom1_prenom1_nom2_prenom2" with all the listed requirements.

The limit date to send your project is: <ins>**4 MAY 2023**</ins>.

Note: It is possible to mix code and report in a notebook, and to mix notebook and independent code. 
However, the report must be contained in a single document.
The criteria for evaluation will be (in order):
- Compliance with instructions,
- Clarity of code,
- Clarity of the report,
- Relevance of the code,
- Results found.

## Last Session:
The last session (12th April 2023) will be a "presentation session":
- 15 minutes
- Present the current state of your work
- Present what you are going to do after (unless you have finished the project)
- Bring your code and draft report.

You don't have to prepare slides or other supports.

Plan to speak ~5 minutes (=> ~10 minutes of Q&A).

## Mandatory general instructions:
- The programming language used must be Python (and only Python),
- No dependencies other than Python modules should be used (Jupiter notebooks accepted),
- Unzipping your archive will produce a folder name "nom_prenom", this is "your folder",
- The code should run once your folder is placed in this repository (i.e. use relative paths only!!),
- The code should not produce any file outside of your folder,
- The code must be commented,
- The code must follow PEP 8 as much as possible,
- The code must contain unit tests,
- The code must be object-oriented.

Note: the objective is to perform an analysis. Therefore, you should not (necessarily) aimed to produce a bioinformatic pipeline robust to changes in the data format for example.
However, it is expected that you ensure that all analyses are correct, in a way or in another.

## Final note:
You will find in the following several steps to help you for the analysis of the data. 
Some steps will tell you exactly what to do, and some will be very open, requiring an investigation and a reflexion. 
Before coding an analysis you could have imagined, present it to me so that I can advise you. 
Also do not rush the avalaible steps, follow the cadence of the sessions. 
For difficult steps, I will initiate and conduct the reflexion. 

# Step 1 - Data Preprocessing
You have access to the raw data of the study in the folder "Data".
This folder contains two information:
- The RNA counts of each sample,
- An annotation file, giving information on the experiment, and (more interestingly) the samples.

Download the data and start to preprocess them.

## Gather RNA counts
To analyze the samples, you will need to merge them into a single Python object.
One standard way to do this is to build a dataframe (or any "table like" struture) such as each row is "sample" and each column is a "gene". Make sure to test your dataset, so that if you change something later on, errors can be catch easily.

Don't forget that the code should be object-oriented.

You can fit your code on the exact state of the data, i.e., design your code to work on thise data on not necessarily on other. 
Here is an example on how you can load the data
```python
import pandas as pd
import glob
import re

path = "./Data" # the path of the data

pdList = [] # variable to temporary store all dataframes (one for each txt file)
# For all txt file
for fname in glob.glob(path+"/*.txt"):
    df = pd.read_table(fname) # put the file in a dataframe
    sample_name = re.search("GSM\d+", fname).group() # search the name of the sample in the file name
    df.rename(index= df["gene/TE"], inplace=True) # rename the index (=rows) using the column containing the gene name
    df.drop(columns=df.columns[0], axis=1, inplace=True) # drop the first column containing the gene name, no more need
    df.rename(columns={ df.columns[0]: sample_name }, inplace = True) # rename the column (there is only one column at this step) using the sample name
    pdList.append(df) # add the current dataframe in the list
data_matrix = pd.concat(pdList, 1) # concat all dataframe in 1 dataframe
data_matrix = data_matrix.transpose() # transpose the dataframe to get a more standard shape (samples x variables)
```
This code should be included in an oriented object code. 
Also, this code assummes that each file is correct, contains no error etc... 
This is something to check to be rigorous. 

## Gather sample annotations
The sample annotations are all placed in a unique "xml" file. First, open the file with any text editor, and try to understand its architecture. Then, identify the information that could be relevant for your analysis. 

Finally, create a dataframe (or any "table like" structure) such as each row is a samples and each column an annotation. Make sur to test your dataset (gene counts+annotations) so that you can catch any error in the next steps.

To parse an xml, you can use the library "xml.etree.ElementTree".
You need to explore the file manually (using a text editor) to catch the structure and the name of all blocks etc...
The samples are contained in blocks named "Sample", and other information are in other blocks that you need to identify.
Here is an example to make a dataframe containing only on column corresponding to the "Cns_subregion".

```python
data_annotation = pd.DataFrame(columns = ['Sample_id', 'Cns_subregion']) # initialisation of the dataframe
xtree = et.parse('./Data/GSE124439_family.xml') # create a variable containing the xml in a tree shape
xroot = xtree.getroot() # get the root of the tree to start the exploration of the tree/xml
# for each element named "sample" that can be found from the root
for child in xroot.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
    temp_sample_id = child.attrib['iid'] # the attribut of this node contains the sample id ()
    # for each element named "Characteristics" that can be found from the current sample
    for child2 in child.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics"):
        if(child2.attrib["tag"] == "cns subregion"):
            temp_cns_subregion = child2.text.replace('\n', '')
    temp_df = pd.DataFrame({'Sample_id': [temp_sample_id], 'Cns_subregion': [temp_cns_subregion]})
    data_annotation = pd.concat([data_annotation, temp_df])
```

## Make first preprocessing functions
By now, you should already have at least one class in your code, with associated getters and setters.
By getters and setters I mean methods that should be used to access to the attributes of instances. 
Good pratice is to prevent the direct access of the instance attributes, so that you can control the way they can be update. 
To do this in python, you can use the prefixe "__" before the attribute name. 
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []
```
This way, the attribute "__data_matrix" can not be change or even read outside of the class. 
A getter (or setter) method need to be developped for this purpose. 
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
```
In this example, the getter method have no real interest, but anyway it is good pratice. 
Think about other preprocessing functions that could usefull for the next steps. 
For example, functions that can check if you have all needed annotations, or function that can subset your data based on some annotation criteria (i.e., get the subdataframe the "control" samples only). Or functions that modified one attribute (e.g., the data matrix) and update accordingly the other attributes (e.g., the annotations).

At last, you may want a convient way to "check" or look at your objects. 
Therefore a "print" function can be defined. 
By simply defining a "__str__" method in your class you can change the behevior of the "print" function when used with your class objects. 
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
    
    def __str__(self):
        return "put here what you want to print"
```

# Step 2 - Descriptive analysis
The descriptive analysis covers all kind of analyses that are direct description of the data, such as computing mean, standard deviation, histogram, boxplot...

## Samples description:
For each sample, compute the mean (across all genes), the median, the standard deviation. Find a way to efficiently report all those data.

Use the annotation data to describe your whole dataset. How many "disease groups", how many sample "sources", how many samples per individual etc... 
This description should guide you for the next step, in order to either correctly group your data when comparing subsets of samples or to avoid potential bias.

You should output all information relative to the samples concisely, so that we have a full (but summarize) view of our sample.

## RNA counts description:
For each gene, compute the mean (across all samples), the median and the standard deviation. Find a way to efficiently report all those data, and make your first interpretation. *Spoiler* : you may have to use graphs and transform/manipulate the data.

Samples correspond to different individuals, to ALS or control individuals etc... Think about what kind of subsets you could analyse and why (and do the descriptive analysis for those subsets).

## Start your report
At this step, you should have already begin your report. 
Before going futher, clean your code and refine your report.
