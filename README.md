# Bioinfo Project 2024

This repository contains materials for the bioinformatics project (2024 L3 INFO, Paris-Saclay). The objective of this project is to identify biomarkers associated with ALS (Amyotrophic Lateral Sclerosis) disease. To achieve this, you will have access to RNA-Seq sequencing data from post-mortem brain cortex biopsies of individuals diagnosed with ALS and those without the disease.

The data for this project originate from the study titled "Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon Activation, Oxidative Stress, and Activated Glia" by Tam et al. The complete study is available [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6866666/).

## Introduction

This README will serve as your guide throughout the project and your analyses. Please follow all the steps listed below. Unless specified otherwise within this document (e.g., mandatory instructions), you are granted complete freedom in choosing how to accomplish each task, particularly in regards to coding techniques.

This document will receive updates after each session. Always refer to the most recent update date mentioned above.

You will work in pairs from the beginning. While discussions with other pairs are encouraged, direct code sharing (such as copy/paste) is strictly prohibited. Throughout the project, meticulously record both your personal contribution and that of your partner. This information will be required during the final session, which will include a roughly 15-minute oral presentation (details to follow).

In the initial sessions, you will receive lessons or introductions on concepts that may be new to you. Once these lessons begin, please pause all other activities and attend the sessions, even if you are already familiar with the concepts being discussed.

At the conclusion of your project, you are required to submit the following components via email to philippe.rinaudo@universite-paris-saclay.fr:
- All your code,
- A comprehensive report detailing the results of your analyses,
- A text file (either .csv, .txt, or similar) containing a ranked list of your top 100 genes (with the most significant gene listed first).

The subject line of the email should read: "Projet_bioinfo_2024_lastname1_firstname1_lastname2_firstname2". Additionally, please include an archive named lastname1_firstname1_lastname2_firstname2" containing all the required elements mentioned above.

The submission deadline for your project is 15 MAY 2024 (to be confirmed).

Note that it is permissible to integrate your code and report within a single notebook, as well as to combine the notebook with standalone code files. However, the report must be consolidated into a single document for submission. The evaluation criteria will prioritize the following, in order:
- Adherence to the provided instructions,
- The clarity of your code,
- The clarity of your report,
- The relevance and efficiency of your code,
- The significance and accuracy of your findings.

## Last Session:

The final session, scheduled for 10th April 2024, will be dedicated to presentations:
- You will have 15 minutes to present,
- Discuss the current progress of your work,
- Outline your next steps, unless your project is already complete,
- Have your code and a draft of your report ready for discussion,
- There is no need to prepare formal slides or other presentation materials,

Aim to speak for approximately 5 minutes, leaving about 10 minutes for questions and answers.

## Mandatory General Instructions:
- Programming Language: Your project must be coded in Python exclusively,
- Dependencies: Only Python modules are allowed for use, including Jupyter notebooks,
- Archive Content: Extracting your archive should create a folder named "nom_prenom" (your folder),
- Code Execution: Ensure your code runs correctly when your folder is placed within this repository by using relative paths only,
- File Management: Your code must not create any files outside of your designated folder,
- Code Comments: All code should be thoroughly commented,
- Coding Standards: Adhere to PEP 8 guidelines as closely as possible for coding style,
- Testing: Include unit tests within your code,
- Design Paradigm: Your code should be designed in an object-oriented manner as much as possible,

Note: The primary goal is to conduct an analysis. You are not required to develop a bioinformatics pipeline that is robust to data format changes. However, it's crucial to ensure the accuracy of all analyses by any necessary means.

## Final Note:
Within this document, you will find a series of steps intended to guide your analysis. Some steps provide clear instructions on what to do, while others are more open-ended, requiring investigation and thoughtful consideration. Before implementing any analyses you conceive, present them for approval to ensure they align with project goals. Do not rush through the available steps; follow the session pace instead. For more challenging steps, I will guide the discussion and thought process.

# Step 1 - Data Preprocessing

You have access to the raw data for the study in the "Data" folder. This folder contains two types of information:
- The RNA counts for each sample,
- An annotation file providing details on the experiment and, more importantly, the samples.
Begin by downloading the data and initiating the preprocessing.

## Gather RNA Counts

To analyze the samples, you will need to consolidate them into a single Python object. A common approach is to create a dataframe (or any "table-like" structure) where each row represents a "sample" and each column represents a "gene". It's crucial to thoroughly test your dataset to ensure that any future modifications can be easily identified and corrected.

Remember, your code must follow an object-oriented design.

Your code can be specifically tailored to the current state of the data, meaning it should be designed to work with this particular dataset rather than being universally applicable. Below is an example of how you can load the data.

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
This code should be integrated into an object-oriented design. Additionally, it assumes that each file is correct and free of errors, which is an aspect that must be verified to ensure thoroughness.

## Gather Sample Annotations

The sample annotations are consolidated into a single "xml" file. Initially, open this file with any text editor to familiarize yourself with its structure. Then, determine which pieces of information are pertinent to your analysis.

Subsequently, construct a dataframe (or any "table-like" structure) where each row represents a sample and each column an annotation. Ensure to rigorously test your dataset (combining gene counts and annotations) to promptly identify and address any discrepancies in subsequent steps.

To parse an XML file, you can employ the "xml.etree.ElementTree" library. It will require manually exploring the file (using a text editor) to understand the structure and identify all relevant sections. The samples are encapsulated within blocks named "Sample," with additional information located in other blocks that you will need to discern.

Below is an example on how to create a dataframe with a single column corresponding to the "Cns_subregion".

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

## Create Initial Preprocessing Functions
At this stage, your code should already include at least one class, accompanied by associated getters and setters. By getters and setters, I refer to methods that facilitate access to an instance's attributes. It is considered best practice to restrict direct access to the instance attributes, allowing you to manage how they can be updated. In Python, this can be achieved by using the prefix "__" before the attribute name, which makes the attribute private and accessible only within its class through defined methods.
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []
```

By adopting this approach, the attribute "__data_matrix" becomes inaccessible for modification or even reading from outside the class. To enable access to this attribute, it is necessary to develop a getter (or setter) method specifically for this purpose. These methods provide a controlled way of accessing and updating the private attribute, ensuring that any changes to "__data_matrix" are handled appropriately within the class's defined interface.

```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
```

In this instance, the getter method might not seem particularly useful, but adopting this approach is considered good practice. It establishes a foundation for more advanced functionality and maintains the integrity and security of the data within your class.

Think about other preprocessing functions that could be beneficial for subsequent steps. For instance, functions that verify the completeness of required annotations or that allow for the subsetting of your data based on specific annotation criteria (e.g., extracting a sub-dataframe for "control" samples only). Additionally, consider functions that can modify one attribute (such as the data matrix) and update related attributes (like the annotations) accordingly.

Lastly, it might be useful to have a convenient method to "check" or examine your objects. Therefore, defining a "print" function can be particularly helpful. This can be achieved by implementing the "str" method within your class, which alters how the "print" function behaves when applied to your class objects. This method allows you to define a custom string representation of your objects, making it easier to understand their state or contents at a glance.

```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
    
    def __str__(self):
        return "put here what you want to print"
```

# Step 2 - Descriptive Analysis

Descriptive analysis encompasses all kinds of direct data descriptions, such as calculating means, standard deviations, and generating histograms and box plots.

## Sample Description:

For each sample, calculate the mean (across all genes), the median, and the standard deviation. Devise an efficient method to report these data comprehensively. Utilize the annotation data to characterize your entire dataset. Detail the number of "disease groups," the variety of sample "sources," and the count of samples per individual, among other aspects. This overview should assist you in planning the next steps, enabling you to group your data correctly when comparing subsets of samples or to identify and mitigate potential biases.

Your goal is to concisely present all relevant information about the samples, providing a comprehensive yet summarized view of your dataset. To achieve this, consider the use of the following visualizations:

- Bar Charts for "Disease Groups" and Sample "Sources": Create bar charts to show the distribution of samples across different disease groups and sources. This will help you quickly visualize the diversity and spread within your dataset.
- Histograms for Measures of Central Tendency and Dispersion: Use histograms to represent the distribution of means, medians, and standard deviations calculated for each sample. This will provide an overview of the variability in gene expression across your study.
- Summary Tables: Prepare tables to summarize key information, such as the number of samples per disease group, source, and individual. Also include a summary table of descriptive statistics (mean, median, standard deviation) for each sample.
- Box Plots: For a more detailed analysis of gene expression distribution, use box plots for each disease group or sample source. This allows for the visualization of medians, quartiles, and outliers.
- Heatmaps: Although a detailed statistical analysis is not required at this stage, a heatmap showing gene expression across different samples or groups can provide a preliminary visualization of expression patterns.

## RNA Counts Description:

For each gene, calculate the mean (across all samples), the median, and the standard deviation. Develop an effective strategy to report these data and offer your initial interpretations. Note that you may need to employ graphical representations and undertake data transformation or manipulation to achieve a clear understanding of the data.

Given that samples represent different individuals, including ALS patients and control subjects, contemplate the subsets of data that could be analyzed and the rationale for selecting these subsets. Conduct the descriptive analysis for these chosen subsets.

To facilitate this analysis, consider the following steps and visualizations:

- Gene Expression Overview: Start with basic statistics (mean, median, standard deviation) for each gene across all samples. This foundational analysis helps identify genes with high variability, potentially indicating differential expression related to disease status or other factors.
- Data Transformation: Given the likely skewed distribution of gene expression data, log transformation (i.e., just log the counts) or other normalization methods might be necessary to stabilize variance and improve the interpretability of statistical analyses.
- Graphical Representations:
-- Box Plots: Display the distribution of expression levels for each gene across ALS patients and control subjects. Box plots are excellent for visualizing the central tendency, dispersion, and outliers within each group.
-- Histograms: Use histograms to illustrate the distribution of expression levels for selected genes, helping to understand the skewness or bimodality of the data.
-- Heatmaps: Create heatmaps to compare gene expression patterns across samples. This is particularly useful for identifying groups of genes that behave similarly across ALS patients versus controls.

## Begin Your Report:

At this juncture, you should have already started drafting your report. Before proceeding further, take this opportunity to refine your code and enhance your report, ensuring clarity and coherence in your presentation of the analyses conducted so far.
