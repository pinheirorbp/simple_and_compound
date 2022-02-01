# Compound topology uncovers complex structure of species interaction networks
Rafael B. P. Pinheiro, Gabriel M. F. Felix and Thomas M. Lewinsohn
<a href="https://zenodo.org/badge/latestdoi/431122160"><img src="https://zenodo.org/badge/431122160.svg" alt="DOI"></a>


In this repository we share code for the entire set of analyses performed in the study. See bellow a brief description of each script file:

(1) topology_test.R : Code for the analyses of the binary networks. The code is prepared for parallel computing and may take a long time to run. This script takes networks in the dataset, calculate network properties, create null models and perform topology tests. The output are often very large files with the raw results.

(2) topology_test_W.R : Code for the analyses of the weighted networks. Same as topology_test, but for weighted networks.

These first two files might be usefull for people interested in applying these same analyses in alternative databases.
Because the outputs of (1) and (2) for our analyses are too large, we do provide them in this repository. They are generated in the folder "/files".

(3) aggregation_results.R: This script aggregate outcomes of (1) and organize the raw data for the analyses in the RMarkdown.

(4) aggregation_results_W.R: This script aggregate outcomes of (2) and organize the raw data for the analyses in the RMarkdown.

Files (3) and (4) generate two "txt" files each in the folder "/results". The first is a table with the results of several analyses (e.g., modularity, nestedness, topology) for each network. The second, whose name always begin with "NETDATA", contains basic informations about the networks (e.g., source, kind of interaction).

(5) binary.Rmd: This file is within the folder "/results" and is used to explore, plot and analyse all binary results. The output is a word document, which is also provided in the folder "/results".

(6) weighted.Rmd: This file is within the folder "/results" and is used to explore, plot and analyse all weighted results. The output is a word document, which is also provided in the folder "/results".

Files (5) and (6) produce a very detailed report of all analyses in our study.

Others:

In the folder "\plots" the scripts "aggregation_results.R" and "aggregation_results_W.R" produce images (.png) of each matrix increasing nestedness visualization and, for modular networks, enhancing the compound topology.<br>
In the folder "\network" all the matrices analyzed in this study are available as txt files.<br>
In the files "IDS_set1.txt", "IDS_set2.txt" and "IDS_set3.txt" we indicate the networks to be analyzed by "topology_test.R". For logistical reasons we analyzed the database in three subsets.<br>
In the files "IDS_set1W.txt" we indicate the networks to be analyzed by "topology_test_W.R".<br>
In the folder "\functions" there are functions that are sourced by the other scripts.<br>
