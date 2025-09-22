# LB2_project_group_4
Laboratory of Bioinformatics 2 - group project \
The members of the group are: Andrea Pusiol, Aurora Mazzoni, Perla Lucaboni and Enrica Cavallo.\\

The steps of the project are: 
* Data collection: retrieve relevant datasets from UniProtKB.
* Data pre-processing: preprocess datasets for cross-validation and benchmarking.
* Analyze and visualize dataset statistics.
* Extract relevant features for classification.
* Implement von Heijne’s algorithm and the SVM classifier.
* Evaluate methods using cross-validation and a blind test set.
* Discuss and report results.
* (Optional) Extend the project with additional predictive methods.
* Prepare a manuscript in the format of a scientific article.



# Data Collection
The [data_collection.ipynb](Data_collection/data_collection.ipynb) script provides TSV and FASTA files for both positive and negative data, containing all preliminary data.\
More in detail, these 4 output files are generated:
* sp_positive.fasta
* sp_positive.tsv
* sp_negative.fasta
* sp_negative.tsv

Retrieved data in detail :

<table>
  <tr>
    <th></th>
    <th>Unfiltered</th>
    <th>Filtered</th>
    <th>TM Helical Proteins</th>
  </tr>
  <tr>
    <th>Positive</th>
    <td>2949</td>
    <td>2932</td>
    <td>/</td>
  </tr>
  <tr>
    <th>Negative</th>
    <td>20615</td>
    <td>20615</td>
    <td>1384</td>
  </tr>
</table>

The results are displayied in the following files :
- [sp_positive.fasta](Data_collection/sp_positive.fasta): FASTA file for the positive dataset 
- [sp_positive.tsv](Data_collection/sp_positive.fasta): .tsv file for the positive dataset 
- [sp_negative.fasta](Data_collection/sp_negative.fasta): FASTA file for the negative dataset 
- [sp_negative.tsv](Data_collection/sp_negative.tsv): .tsv file for the negative dataset 
