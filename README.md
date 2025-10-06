# LB2_project_group_4

The members of the group are: Andrea Pusiol, Aurora Mazzoni, Perla Lucaboni and Enrica Cavallo.

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
Our objective is to collect datasets that can be used to train predictive methods for identifying the presence of signal peptides in proteins.
For this project, we will restrict our analysis to eukaryotic proteins, deliberately excluding bacterial and archaeal sequences.
To this end, we require two types of datasets:
* **Positive dataset**: eukaryotic protein sequences that contain experimentally validated signal peptides.
* **Negative dataset**: eukaryotic protein sequences that do not contain a signal peptide.
In both cases, the presence or absence of a signal peptide is annotated in UniProt under the feature “PTM-Processing / Molecule Processing / Signal peptide”.

To build a preliminary positive dataset, we will apply the following selection criteria in the UniProt Advanced Search:
* exclude protein fragments `fragment:false`
*	include only eukaryotic proteins `taxonomy_id:2759`
* select proteins longer than 40 amino acids `length:[40 TO ]`
*	restrict to reviewed (Swiss-Prot) entries `reviewed:true`
*	retain proteins with experimentally validated signal peptides `ft_signal_exp:`

Additionally, in our script [data_collection.ipynb](Data_collection/data_collection.ipynb) we will further refine the dataset by keeping only:
* proteins with signal peptides longer than 14 amino acids
*	entries with protein-level evidence
*	proteins with an annotated cleavage site

To construct a preliminary negative dataset, we will apply the following selection criteria in UniProt:
*	exclude protein fragments `fragment:false`
*	include only eukaryotic proteins `taxonomy_id:2759`
*	select proteins longer than 40 amino acids `length:[40 TO ]`
*	exclude any protein annotated with a signal peptide (at any evidence level) `ft_signal:`
*	restrict to reviewed (Swiss-Prot) entries `reviewed:true`
*	retain only proteins with protein-level existence evidence `existence:1`
*	include only proteins experimentally verified to be localized in one of the following compartments: cytosol, nucleus, mitochondrion, plastid, peroxisome, or cell membrane 
`((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))`

After the selection step, both the positive and the negative datasets must be further filtered and then stored in two different formats.

The first format is a TSV file, which summarizes the most relevant information about the proteins included in each dataset. The reported fields differ between positive and negative datasets.

For the positive dataset, the TSV file [sp_positive.tsv](Data_collection/sp_positive.fasta) will include:
-	UniProt accession number of the protein
-	Organism name
-	Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
-	Protein length
-	Position of the signal peptide cleavage site

For the negative dataset, the TSV file [sp_negative.tsv](Data_collection/sp_negative.tsv) will include:
-	UniProt accession number of the protein
-	Organism name
-	Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
-	Protein length
-	Presence of a transmembrane helix within the first 90 residues (reported as true or false).
  Information about transmembrane helices is particularly relevant, since signal peptides are physico-chemically similar to transmembrane segments. For this reason, we
  expect a higher false positive rate in transmembrane proteins, and including this feature will be useful for evaluating prediction performance.
 	
The second format is a FASTA file which reports the protein sequences.
- [sp_positive.fasta](Data_collection/sp_positive.fasta): FASTA file for the positive dataset
- [sp_negative.fasta](Data_collection/sp_negative.fasta): FASTA file for the negative dataset 


In summary :

|             | Unfiltered | Filtered | TM Helical Proteins |
|:-----------:|:----------:|:--------:|:-------------------:|
| **Positive**|   2949     |   2932   |         /           |
| **Negative**|  20615     |  20615   |        1384         |






# Data Preparation

Once the two preliminary datasets, positive and negative, have been obtained, they must be pre-processed in parallel.
The first step is to remove redundant sequences from each dataset, in order to obtain non-redundant datasets.
This procedure is essential for several reasons:
* **Reduction of redundancy**: it prevents protein families or organisms that are overrepresented, as often happens in UniProt databases, from excessively biasing the model.
* **Balanced dataset**: it ensures that each family, fold, or protein function has a more balanced weight, allowing for a more accurate and fair analysis.
* **Prevention of data leakage**: by removing highly similar sequences, it avoids the risk that nearly identical copies end up in both the training and benchmarking sets, a situation that would lead to artificially inflated performance and an unrealistic evaluation of the model’s generalization ability.

This is done through the process of sequence clustering, which groups sets of similar proteins using two fundamental parameters: sequence identity and alignment coverage. For two sequences to be considered part of the same cluster, both conditions must be satisfied.

To perform this process, we use MMSeqs2, one of the most efficient tools for large-scale clustering. The program compares all sequences, groups them into clusters according to the chosen parameters, and selects a single representative for each cluster, thereby eliminating redundant sequences. The command to run sequence clustering with MMSeqs2 follows the general syntax:

```sh
mmseqs easy-cluster input.fa cluster-results tmp --min-seq-id 0.3 \ -c 0.4 --cov-mode 0 --cluster-mode 1
```

Where:

The `input.fa file` represents the initial dataset in FASTA format.

The clustering results are saved with the prefix `cluster-results`.

`tmp` is a temporary folder automatically created to store intermediate files during execution.

`--min-seq-id 0.3` sets a minimum threshold of 30% sequence identity: this means that two sequences can be considered similar only if at least 30% of their residues are identical.

`-c 0.4` imposes a coverage threshold of 40%, meaning that at least 40% of the sequence length must actually be aligned.

`--cov-mode 0` specifies that the coverage requirement must be satisfied for both the query sequence and the target sequence, which is to say in both directions.

`--cluster-mode 1` selects the clustering algorithm, in this case the connected component mode, which builds clusters as sets of sequences connected to each other within a similarity network.

When we perform sequence clustering with MMSeqs2, we mainly obtain two output files:

* The first is **cluster-results_rep_seq.fasta**, a FASTA file that contains all the representative sequences, one for each cluster identified. ([cluster-results_positive_rep_seq.fasta](Data_preparation/cluster-results_positive_rep_seq.fasta), [cluster-results_negative_rep_seq.fasta](Data_preparation/cluster-results_negative_rep_seq.fasta)) 

* The second file is **cluster-results_cluster.tsv**, in tabular format. It contains two columns: in the first column, the ID of each original sequence from the input file is reported, while the second column shows the ID of the representative sequence of the cluster to which that sequence has been assigned. In this way, we obtain a kind of “map” that tells us, for each input sequence, which cluster it belongs to and what its representative is. ([cluster-results_positive_cluster.tsv](Data_preparation/cluster-results_positive_cluster.tsv), [cluster-results_negative_cluster.tsv](Data_preparation/cluster-results_negative_cluster.tsv))

Once redundancy has been removed, each dataset is divided into two distinct subsets:

* **Training set** (about 80% of the data), used to train the models, optimize hyperparameters, and perform cross-validation experiments.

* **Benchmarking set** (about 20% of the data), used exclusively in the final step to evaluate the model’s generalization ability.

Finally, the training set of each reduced dataset is divided into 5 subsets.
The split was performed randomly, while ensuring that each subset preserved the same ratio of positive and negative sequences as in the original dataset.
This allows the implementation of 5-fold cross-validation, a procedure that enables training and validating the model multiple times, by iteratively using 4 folds for training and 1 fold for validation.
To ensure traceability, the fold assignment of each protein was also recorded, so that the data partitioning can be clearly and transparently reconstructed.



Output: [Positive](Data_preparation/positive_set.tsv) and [Negative](Data_preparation/negative_set.tsv)

In summary:

| **Dataset**   | Positive | Negative |
|:-----------------------:|:--------:|:--------:|
| **Before Clustering**   |   2932   |  20615   |
| **Cluster representatives** |  1093   |   8934   |
| **Training**            |     875     |     7147     |
| **Test**                |     218     |     1787     |






# Data Analysis

The main operations of this phase are:

- importing the positive and negative datasets;

- adding labels to distinguish positive and negative sequences;

- checking the uniqueness of protein identifiers;;

- setting `protein_id` as the index of the tables;

- renaming the  `class ` column to  `dataset_class ` to distinguish training and test;

- concatenating the two datasets into a single dataframe;
  
- final separation into training and test datasets based on the `dataset_class` column.

Once the dataset has been structured, before using it to train predictive models it is essential to carry out a preliminary statistical analysis. The purpose of this step is to verify that the available data are of adequate quality, balanced, and truly representative of the biological phenomenon we are going to study.

The first step consists of computing basic statistics on the training and benchmark datasets. This includes measuring values such as mean, median, standard deviation, protein length distributions, or amino acid frequencies. These indicators provide a quantitative overview of the dataset and help identify potential anomalies, imbalances, or outliers.

Next, it is necessary to assess the adequacy of the dataset with respect to the goals of the analysis. One must ask whether the number of sequences is sufficient, whether the positive and negative classes are well balanced, and whether the categories of interest (organisms, functions, lengths) are properly represented. This verification is crucial to ensure that the trained models can produce reliable results that are not biased by partitioning artifacts.

Finally, it is important to demonstrate that the datasets represent a true “snapshot” of the real protein space. In other words, the data should reflect the variability and the typical features of natural proteins, without being limited or distorted by the way they were extracted. Only datasets that faithfully capture biological diversity can lead to generalizable and biologically meaningful results.

To achieve these objectives, several descriptive plots are generated:

- We examined the **distribution of protein lengths** , comparing the positive dataset (proteins with signal peptides, secretory pathway) and the negative dataset (proteins without signal peptides). The data were visualized as normalized histograms with overlaid density curves, shown separately for the training set and the test set.
The goal of this step was twofold. On one side, it allows us to assess whether protein length could serve as a discriminative feature between secretory and non-secretory proteins. If consistent differences are observed between positive and negative distributions, this property might contribute useful information to classification models. On the other hand, if the two distributions largely overlap, protein length is unlikely to represent a strong predictive feature.
By plotting the training and test sets separately, we also ensure that the data split did not introduce any bias. The distributions appear consistent across the two subsets, confirming that the benchmark set is representative of the same underlying data space as the training set.






# VonHeijne method


| **Metric**   |      **Mean ± SE**       |
|------------- |----------------------|
| Accuracy     | 0.936  ± 0.0032      |
| Precision    | 0.7671 ± 0.0311      |
| Recall       | 0.6023 ± 0.0198      |
| F1 Score     | 0.6743 ± 0.0132      |
| MCC          | 0.6454 ± 0.0155      |
| Threshold    | 9.2161 ± 2.6279      |











