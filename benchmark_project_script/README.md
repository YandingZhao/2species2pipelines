# benchmark_cross_species_integration

In the study, we benchmark nine methods (fastMNN, harmony, BBKNN, SATURN, SAMap, seurat v4 CCA, scGen, scVI, scanorama) for cross-species integration 

The relevant files should be provided in the following structure in order to running the code successfully.
If not, then the user need to manuelly change the saving/reading path in the code

note: the folders for each method under method_output and evaluation should be named as below separately:
fastmnn, bbknn, harmony, saturn, samap, scanorama, scGen, scVI, seurat4

|----script # folder for the script, which is where we are now
|
|          |-method_outputs # store the integration output for each method separately (subfolers for each method)
|----output|               
|          |-evaluation  # store the evaluation result for each method separetely (subfolders for each method)
|          |
|          |-final_evaluation_table # final table for each methods' all evaluation results 
|
|----data/process_data|-h5ad_files # this is to store the preprocessed h5ad file
|                     |-rds_files # this is to store the preprocessed rds file 
|
|----OrthoFinder/one2one_orthologs  # this is the folder store the one2one mapping table for each task

We followed the following step for the evaluation:
Step 0: environment install
we provided three environment for the benchmark and evaluation in the conda_env folder
specifically:
benchmark_R.yml is for the integration methods using R language
benchmark_py.yml is for the integration methods using python language
eval_scib.yml is for the evaluation using scib package

Step 1: method integration
In this step, we perform the integration based on each method.
Since Seurat v4, Harmony, and fastMNN are executed in R, while the other methods are executed in Python, we used different environments to avoid package incompatibilities. Beside, in the top of each code, we provide comment on what the code is doing and how to run the code.

For the methods run in R, please use the benchmark_R environment (benchmark_R.yaml)
    under this environment, we need to run the following script for R-based methods integration
    seurat4_for_all_tasks.R
    harmony_for_all_tasks.R
    fastMNN_for_all_tasks.R

    usage for all these scripts are similar,we just need to specify the task name to the script.
    For example: Rscript seurat_for_all_tasks.R task1 task2 

For the methods (bbknn, scVI, scGen, scanorama) run in python, please use the benchmark_py environment (benchmark_py.yaml)
    under this environment, we need to run the following scripts for python-based methods integration
    bbknn_for_all_tasks.py
    scVI_for_all_tasks.py
    scGen_for_all_tasks.py
    scanorama_for_all_tasks.py

    usage for all these scripts are similar,we just need to specify the task name to the script.
    For example: python scVI_for_all_tasks.py task1 task2 

For the methods (saturn, samap) have complicated running processes in python, please follow the official tutorial.
    for saturn, please refer: https://github.com/snap-stanford/SATURN/blob/main/Vignettes/frog_zebrafish_embryogenesis/Train%20SATURN.ipynb
    for samap, please refer: https://github.com/atarashansky/SAMap/blob/main/SAMap_vignette.ipynb

Step 2: Evaluation on method outputs
In this step, we evaluate the methods' output in 13 metrics: ARI(batch), ARI(celltype), ASW(batch), ASW(celltype), NMI(batch), NMI(celltype), iLISI, cLISI, kBET, PCR batch, graph connectivity, Highly variable gene conservation, trajectory conservation.

specifically, for the ARI, ASW, NMI and LISI evaluation on the non-graph output methods (which means except samap, bbknn), please run the following script:
    emb_recon_ARI_LISI.R (benchmark_R environment)
    emb_recon_ASW_NMI.py (benchmark_py environment)

for the bbknn and samap evaluation, please run the below script in order:
    (1) generate_leiden_clusters.py (benchmark_py environment)
    (2) samap_bbknn_ARI.r (benchmark_R environment)
    (3) samap_bbknn_evaluate_NMI.py (benchmark_py environment)
    (4) samap_bbknn_evaluate_LISI.py (benchmark_py environment)

for the kBET, PCR batch, graph connectivity, Highly variable gene conservation, trajectory conservation evaluation on all methods, please run below:
    scib_metric_running.py (eval_scib environment)

To gather all the evaluation results, we can run the following to generate the final table for visualization.
    generate_final_evaluation_result.py (benchmark_py environment)

Step3: visualization for the final evaluation results
please run the code below in R
    visualization.R

Step4 (optional): for cell type tree construction based on integrated output, please follow the official tutorial:
https://github.com/dunnlab/cellphylo/tree/master/analysis

### extra note about how to run orthofinder to get the one2one mapping table
In our study, we used OrthoFinder to identify one-to-one orthologous genes between species. OrthoFinder is a robust tool designed for comprehensive orthology analysis across multiple species. It identifies orthogroups, which are sets of genes descended from a single gene in the last common ancestor of all the species under consideration. To run orthofinder, please follow the official tutorial: https://github.com/davidemms/OrthoFinder


Specifically, OrthoFinder determines one-to-one orthologs through a series of steps: 
OrthoFinder determines one-to-one orthologs through a series of steps: 
(1) All-versus-All Sequence Similarity Searches: 
a. It performs all-versus-all BLAST searches to calculate sequence similarities between genes from all species. 
b. These similarities are normalized to account for gene length and evolutionary distance, reducing biases in the similarity scores. 

(2) Clustering into Orthogroups: 
a. Using the normalized similarity scores, genes are clustered into orthogroups using the Markov Cluster Algorithm (MCL). 
b. Each orthogroup contains genes that are all descended from a single ancestral gene. 

(3) Gene Tree Construction: 
a. OrthoFinder constructs gene trees for each orthogroup. 
b. These trees represent the evolutionary relationships of the genes within an orthogroup. 

(4) Species Tree Inference: 
a. It infers an unrooted species tree based on the gene trees using STAG algorithm. 
b. The species tree is then rooted using STRIDE to identify the direction of evolution. 

(5) Ortholog Identification: 
a. By comparing gene trees to the rooted species tree, OrthoFinder identifies orthologous relationships. 
b. One-to-one orthologs are identified where a single gene in one species is directly descended from a single gene in another species without any duplication events. 

For more detail information about the method, please read their original paper:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y

