# LB2_project_group_4
Laboratory of Bioinformatics 2 - group project \
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


Retrieved data in detail :

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
`--min-seq-id 0.3` sets a minimum threshold of 30% sequence identity: this means that two sequences can be considered similar only if at least 30% of their residues are identical.
`-c 0.4` imposes a coverage threshold of 40%, meaning that at least 40% of the sequence length must actually be aligned.
`--cov-mode 0` specifies that the coverage requirement must be satisfied for both the query sequence and the target sequence, i.e., in both directions.
`--cluster-mode 1` selects the clustering algorithm, in this case the connected component mode, which builds clusters as sets of sequences connected to each other within a similarity network.


The output consists in 2 files : \
Cluster-results_rep_seq.fasta → a FASTA file containing all the representative sequences, one for each found cluster. ([positive](Data_preparation/cluster-results_positive_rep_seq.fasta), [negative](Data_preparation/cluster-results_negative_rep_seq.fasta)) \
Cluster-results_cluster.tsv → a TSV containing two columns: reports the ID id each sequence in the input file & reports the ID of the representative sequence identifying the cluster the sequence in column 1 has been assigned to. ([positive](cluster-results_positive_cluster.tsv), [negative](Data_preparation/cluster-results_negative_cluster.tsv))


Now the .fasta files containing only "unrelated" proteins will be used to split each dataset into two separate subsets:

Training set: used to train the methods, optimize model hyperparameters and perform cross-validation experiments. \
Test set: used to test the generalization performance of the different models. 

For this purpose we developed the data_preparation.ipynb file, that takes as input:
Cluster-results_rep_seq.fasta files
sp_positive/negative.tsv files

Output: [Positive](positive_set.tsv) and [Negative](negative_set.tsv)


| **Dataset**   | Positive | Negative |
|:-----------------------:|:--------:|:--------:|
| **Before Clustering**   |   2932   |  20615   |
| **Cluster representatives** |  1093   |   8934   |
| **Training**            |     875     |     7147     |
| **Test**                |     218     |     1787     |


# Data Description

At this point, we need to perform statistical analysis of the datasets, including distributions of different aspects of data, like compositions, taxonomy and SP lenght. This is usefull to avoid biases in the data. Furthermore, the results of these analysis will be represented using different types of plots, in order to describe the training and the benchmarking datasets indipendently.
Specifically, we will produce the following graphs:

* The distribution of protein lenghts, comparing positive and negative sequences

* The distribution of SP lenghts

* Comparative amino-acid composition of SPs against some background distribution

* Taxonomic classification, at kingdom and species level

* Sequence logos of SP cleavage sites, extracting the cleavage-site motifs [-13,+2]. This is a typical representation of sequence conservation, starting from aligned sequences. The logo consist of a stack of letters for each aligned position, where the height of the entire stack informs you about the information content, while the height of each letter is proportional to its frequency at that position.



# VonHeijne method




| **Metric**   | **Value** |
|--------------|-----------|
| Accuracy     |           |
| Precision    |           |
| Recall       |           |
| F1 Score     |           |
| MCC          |           |










