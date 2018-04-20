# &#8734;SCITE
========

## Overview
--------------


**&#8734;SCITE** is a variant of the software package SCITE for the computation of mutation histories of somatic cells. It has two new main features compared to the original software:

* **Doublets**: &#8734;SCITE comprises a model for doublets, i.e. samples with mutation profiles that come from sequencing two cells instead of just one.
* **Testing the ISA**: &#8734;SCITE can be used to test the validity of the infinite sites assumption for a single-cell mutation matrix. 



Given the noisy mutation profiles of single cells, **&#8734;SCITE** (like SCITE) performs a stochastic
search to find the Maximum Likelihood (ML) or Maximum aposterori (MAP) tree and/or to sample from the posterior
probability distribution. Tree reconstruction can be combined with an estimation
of the error and doublet rates in the mutation profiles.

**&#8734;SCITE** allows the user to **specify a single mutation as recurrent**. It then searches a tree space that represents M by two separate nodes. (The second node represents either an independent acquisition of the mutation in another branch, or the loss of the mutation if the two nodes occur in the same lineage.)

To perform **model selection** between the infinite sites model (allowing no recurrences or losses) and the finite sites model allowing recurrences you will need the run **&#8734;SCITE** with each possible (single) recurrence and without any recurrence.
More details on the doublet model and the model selection can be found in:

###### Citation:

_Kuipers, J., Jahn, K., Raphael, B. J., & Beerenwinkel, N. (2017). Single-cell sequencing data reveal widespread recurrence and loss of mutational hits in the life histories of tumors. Genome research, 27(11), 1885-1894._

## Availability
---------------

**&#8734;SCITE** is freely available under a GPL3 license at https://gitlab.com/jahnka/infiniteSCITE

##    How to run **&#8734;SCITE**
--------------------------



### Mac OS X

To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type

	clang++ *.cpp -o infSCITE

This writes a file named `infSCITE`. With some compiler versions you may need to use the option `-std=c++11`.

Assuming the sample data file pat_2.csvis located in the same folder, `infSCITE` can then be executed as follows

	./infSCITE -i pat_2.csv -n 16 -m 115 -r 1 -l 700000 -fd 3.45e-3 -ad 1.46e-1 -s -e .2 -p 10000 -d -rec 13

This call returns the MAP tree and samples from the posterior distribution for the given dataset under the assumption that the mutation at position 13 has a recurrence. With this call **&#8734;SCITE** also learns the false negative and doublet rates from the data. See below for other program options.

### Linux/Unix

To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type

	g++ *.cpp -o infSCITE
	
This writes a file named `infSCITE`. With older compiler versions you may need to use the option `-std=c++11`.

Assuming the sample data file pat_2.csvis located in the same folder, `infSCITE` can then be executed as follows

	./infSCITE -i pat_2.csv -n 16 -m 115 -r 1 -l 700000 -fd 3.45e-3 -ad 1.46e-1 -s -e .2 -p 10000 -d -rec 13

This call returns the MAP tree and samples from the posterior distribution for the given dataset under the assumption that the mutation at position 13 has a recurrence. With this call **&#8734;SCITE** also learns the false negative and doublet rates from the data. See below for other program options.


##  Input Files
---------------


### 1. Mutation Matrix


Each column specifies the mutation profile of a single cell, and each row
represents one mutation. A white space character separates the columns.

#### (a) Only absence/presence of mutation is distinguished
The entry at position [i,j] should be

* 0 if mutation i is not observed in cell j,
* 1 if mutation i is observed in cell j, or
* 3 if the data point is missing

The sample datasets dataNavin.csv and dataXu.csv have this format.	
#### (b) Heterozygous and homozygous mutations distinguished
The entry at position [i,j] should be

* 0 if mutation i is not observed in cell j,
* 1 if heterozygous mutation i is observed in cell j
* 2 if homozygous mutation i is observed in cell j
* 3 if the data point is missing

The sample datasets dataHou18.csv and dataHou78.csv have this format

### 2. Mutation names (optional)


A list specifying the names of the mutations, e.g. the name of the gene in which
the mutation occurs. For the sample datasets provided here, these files have the extension *.geneNames*. If no such file is specified, the mutations are numbered from 1 to n.

### 3. The true tree (optional)


If the true mutation tree is known (as in simulation experiments), the file with
the true tree (in GraphViz format) can be specified to compare the
predicted trees internally with the true tree.

##  Output Files
----------------

### 1. ML/MAP trees

ML/MAP trees are written to files in GraphViz and Newick format. Files are numbered
consecutively (e.g. dataHou18_ml1.gv, dataHou18_ml1.newick, dataHou18_ml2.gv, dataHou18_ml2.newick, ...). The base name of the output file is derived from the name of the input file (unless a different name is specified via `-o <filename>`).

### 2. Samples from the posterior distribution (optional)

When the `-p <INT>` option is set, **&#8734;SCITE** samples from the posterior distribution, and writes the sampled trees (in parent vector format) together with their scores and learned error rates to a single file (one sample per line). The name of the output file is derived from the input file name using the ending *.sample*.



## Parameters
-------------

*	`-i <filename>`     Replace \<filename\> with the file containing the mutation matrix

*	`-n <INT>`  Replace \<INT\> with the number of mutations (rows) in the dataset.

* `-m <INT>`  Replace \<INT\> with the  number of cells (columns) in the dataset.

* 	`-r <INT>`  Set \<INT\> to the desired number of repetitions of the MCMC.

* 	`-l <INT>`  Set \<INT\> to the desired chain length of each MCMC repetition



##### In case only absence/presence of a mutation is distinguished

*	`-fd <DOUBLE>` Set \<DOUBLE\> to the estimated false positive rate (false discoveries) of the sequencing experiment.

*	`-ad <DOUBLE>` Set \<DOUBLE\> to the estimated false negative rate (allelic dropout) of the sequencing experiment.

##### In case heterozygous and homozygous mutations are distinguished

*	`-fd <DOUBLE>` Set \<DOUBLE\> to the estimated false positive rate (false calling of heterozygous mutation) of the sequencing. experiment

*	`-ad <DOUBLE> <DOUBLE>` Set the first \<DOUBLE\> to the estimated rate of missed heterozygous mutations in the sequencing experiment and set the second \<DOUBLE\> to the estimated rate of heterozygous mutations called as homozygous mutations (dropout of the normal allele).

*	`-cc <DOUBLE>` Set \<DOUBLE\> to the estimated rate of non-mutated sites called as homozygous mutations.   



## Optional parameters
----------------------

#### Program variants


*	`-s` Setting this option causes the sample attachment points (i. e. where the cells would attach to the tree) to be marginalized out.

*	`-p <INT>` When setting this option, **&#8734;SCITE** samples from the posterior distribution, and writes the trees to a file using the parent vector format. The value of \<INT\> specifies how dense the sampling is. The name of the output file is derived from the input file name using the ending *.sample*.
To make sure that **&#8734;SCITE** samples from the posterior distribution `-p <INT>` needs to be combined with `-s` (single cell attachment to the tree is marginalized out) and `-g 1`, or without the `-g` option (gamma is set to 1 by default).

	__Sample file columns (no error rate learning):__

	*	tree log score
	*	number of branches in tree
	*	parent vector

	__Sample file columns (no error rate learning):__

	*	tree log score
	*	number of branches in tree
	*	beta
	*	alpha
	*	combined log score (tree and theta)
	*	relevant doublet rate
	*	doublet rate
	*	parent vector  

*	`-transpose`	This changes the tree representation from mutation tree to rooted binary leaf-labeled tree, where the samples are the leaf labels and mutations are placed on the edges. Using this option can decrease the search space size when there are more mutations than samples. This option is incompatible with the `-s` mode and not implemented for doublets and recurrent mutations.

#### Doublet model

*	`-d`	Invokes the use of the doublet model. Each sample will be treated as a mixture of a singlet and a doublet sample. If no doublet rate is specified, it is optimized for each tree. 

*	`-d <DOUBLE>`	It is possible to use a fixed doublet rate, by specifying it after the `-d` option. No per tree optimization is performed in this case.

Note: Using the `-d` or `-d <DOUBLE>` option makes the program **significantly slower**.
This option is only available for the mutation tree model (default option), not the binary tree model.


#### Recurrent mutation

* `-rec <INT>` Replace \<INT\> with the number of the mutation that should be represented twice in the tree. The number of a mutation is based on the order of occurrence in the mutation matrix. When used for testing the infinite sites assumption, we recommend to always combine this option with the doublet model option `-d`.
 

#### Error learning

*	`-e <double>`   Invokes the learning of error rates. Set \<double\> to a value between zero and one to specify the probability of a move for changing the error rate in the MCMC.

*	`-x <double>`   Scaling of the known error rate for the MH jump (default is 10).


###### Learning beta
When option `-e <double>` is used by default only the false negative rate is learned. (See below for learning alpha.)

Optional parameter:

* `-sd_beta <DOUBLE>`	Sets standard deviation of prior for false negative error rate (default is 0.1). For the mean, the error rate specified after `-ad` is used.

###### Learning alpha

To learn the false positive rate alpha, it is necessary to take into account all sequenced positions not just the ones used in the tree inference. For example when all sequenced positions with normal state in all cells and positions only mutated in one cell are discarded, we have the following additional 0s and 1s:

* z0 = #zeroes in the discarded part of the matrix. In the above example this would be: 
*m \* \#number of all zero columns + (m-1) \* \#number of columns with a single one*

* z1 = #ones in the discarded part of the matrix. In the above example this would be: *#number of columns with a single one*

To invoke the learning of alpha, set:

* `-z <INT> <INT>` Replace the first \<INT\> by z0 and the second \<INT\> by z1

Note: Setting `-z 0 0` is possible but will lead to an overestimation of alpha. If you don't know how to set z0 and z1 but want to still learn alpha, do not pass any values:

* `-z`  (This setting will use some default values for z0 and z1.)

Optional parameter:

* `-sd_alpha <DOUBLE>`	Sets standard deviation of prior for false positive error rate (default is 0.1). For the mean, the error rate specified after `-fd` is used.


#### Output

*	`-o <filename>`   Replace \<filename\> with the desired base of the output file to overwrite the default output file names.

*	`-names <filename>` Replace \<filename\> with a file listing the mutation names. By default the mutations are numbered from 1 to n by order of appearance in the inputfile.

*	`-a` When setting this option, **&#8734;SCITE** adds the individual cells as additional nodes (leafs) to the reported trees. Cells are attached where they fit best in the given tree (with respect to the error rates). By default, only the mutation tree is reported.

*	`-max_treelist_size <INT>`	 This limits the number of co-optimal trees written to output files to the specified <INT> value.

#### Other parameters

*	`-g <DOUBLE>` For ML/MAP computation only: Set \<DOUBLE\> to the desired value of gamma (gamma > 1: more local exploration, possibly local optimum; gamma < 1: easier to explore the space, but less deeply). The default value of gamma is 1 which is necessary for the MCMC chain to converge to the posterior distribution.

*	`-seed <INT>`   Replace \<INT\> with a positive integer to be used as a fixed seed for the random number generator.

*	`-no_tree_list`		This turns off the collection of optimal trees. It can be used to speed up the program when only sampling from the posterior distribution.

*	`-t <filename>`  Replace \<filename\> with a file containing the true tree in GraphViz format.

*	`-move_probs <double> <double> <double>`   Changes the default probabilities for the three MCMC moves to the spedified values. The first move is *prune and re-attach*, the second is *swap node labels*, the third is *swap subtrees*. The default values are (0.55, 0.4, 0.05).

When combined with `-transpose` there are only two move types, *prune and re-attach* and *swap leaf labels* with default probabilities (0.4, 0.6).  Therefore the parameter format for changing these move probabilities is `-move_probs <double> <double>`.


#### Source files

* findBestTrees.cpp
* doublets.cpp
* doublets.h
* recMut.cpp
* recMut.h
* mcmc.cpp
* mcmc.h
* mcmcBinTreeMove.cpp
* mcmcBinTreeMove.h
* mcmcTreeMove.cpp
* mcmcTreeMove.h
* scoreBinTree.cpp
* scoreBinTree.h
* scoreTree.cpp
* scoreTree.h
* trees.cpp
* trees.h
* treelist.cpp
* treelist.h
* rand.cpp
* rand.h
* matrices.cpp
* matrices.h
* output.cpp
* output.h


* findBestTrees.cpp doublets.cpp recMut.cpp mcmc.cpp mcmcBinTreeMove.cpp mcmcTreeMove.cpp scoreBinTree.cpp scoreTree.cpp trees.cpp treelist.cpp rand.cpp matrices.cpp

-i "/Users/jahnka/Desktop/MykolaProject/data/TNBC/TNBC.txt" -n 547 -m 16 -r 2 -l 500 -g 1 -fd 6.04e-6 -ad 0.11545 0.11545 -cc 1.299164e-05 -names "/Users/jahnka/Desktop/MykolaProject/data/TNBC/TNBC.geneNames" -e 0.1 -transpose


g++ findBestTrees.cpp doublets.cpp recMut.cpp mcmc.cpp mcmcBinTreeMove.cpp mcmcTreeMove.cpp scoreBinTree.cpp scoreTree.cpp trees.cpp treelist.cpp rand.cpp matrices.cpp output.cpp binTree_output.cpp -std=c++11 -o SCITE








