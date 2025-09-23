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

# Data Preparation

<table>
  <tr>
    <th></th>
    <th>Total</th>
    <th>Cluster representative</th>
  </tr>
  <tr>
    <th>Positive</th>
    <td>2932</td>
    <td>1092</td>
  </tr>
  <tr>
    <th>Negative</th>
    <td>20615</td>
    <td>8933</td>
  </tr>
</table>






