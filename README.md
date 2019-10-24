# CDMSuite
Community Dectection algorithms applied to PPI network data

------------------------------------------------------
### Community Detection using Modulairty (2012)    ###
### Suite of Modualrity based community detection  ###
### algorithms                                     ###
### Dr. Colin Mclean <cmclean5@staffmail.ed.ac.uk> ###
------------------------------------------------------

Based on the fast community algorithm implemented by Aaron Clauset <aaron@cs.unm.edu>, Chris Moore, Mark Newman, and the R IGraph library Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>.

Suite containing three community detection algorithms based on the Modularity measure containing: Geodesic and Random Walk edge Betweenness [1] and Spectral Modularity [2]. 

----------------------------
*** TO INSTALL AND BUILD ***
----------------------------

(1) Build clustering package within a convienient location,  by running the Makefile, we compiled using gcc version 4.4.7. This should create the executable 'run'. 

(2) To use OpenMP run MakeFile with USEOMP option set to 1:

 > make clean
 > make USEOMP=1 (    to use OpenMP)
 > make USEOMP=  (not to use OpenMP) 

----------------------------
***       TO RUN         ***
----------------------------

To run the clustering package at the command line type the following:

 > ./run  -file testData/karate.txt -a 1 -cols 3

This will execute the Geodeic edge Betweenness algorithm on the network file. The -file argument specifies the network file to analysised. The -a argument selects which of the three algorithms to execute: 1 = Geodesic egde Betweenness, 2 = Random Walk edge Betweenness, 3 = Spectral Modularity. The -cols argument indicates the number of columns in the network file, i.e. whether the network file is weighted (3) or unweighted (2).

The network file can have single or multiple header lines at the begining, which can be skipped using the -header or -skip arguments. The first two tab separated columns in the network file can contain the alpha-numerical IDs for the network edges, for example gene names or Entrez gene IDs. An additional third column can contains the edge weights. Output from running the clustering algorithm can be found in the folder OUT.

 # Network file Example 1
 
 10006 10458
 
 7316 10152
 
 6709 10006

 # Network file Example 2
 
 ADA	POTEF
 
 ADA	ACTB

 BEND7	TRIM27
 
 BEND7	C1QBP
 
 BEND7	LRP4
 
 CCNI2	HSPA1B
 
 CCNI2	HSPA8

The use can obtain help on running each of the arguments by typing the following:
 
 ./run -help 

 ./run requires at least 3 arguments:
 
 -file                     : the network file to run
 
 -header                   : use if there's a file header
 
 -skip N <int>             : use to skip N lines from start of network file
 
 -quite                    : turn off comments to command line [short cut -q]
 
 -seed N <int>             : Set N for the random number generator seed; default is 1.
 
 -a N (1,2,3) <int>        : Set N for the type of algorithm to run; default is 1. Where:
 
                           : 1 = Geodesic edge Betweenness
                           
                           : 2 = Random edge Betweenness
                           
                           : 3 = Spectral Betweenness
                           
 -weighted                 : Specify if network file is weighted (network file with 3 columns)
 
                           : or not (network file with 2 columns ) [sort-cut -w]
                           
 -cols N (2,3) <int>       : Specify if network file is weighted (N=3) or unweighted (N=2)"
 
 -subsample                : Subsample the node set, ramdonly selecting -per N [0.0,1.0] of the node.
 
 -per N [0.0,1.0] <double> : Randomly selected N % of the node set.
 
 ------------------------------------------------------------------
 Example      : ./run  -file testData/karate.txt -a 3 -cols 3
 ------------------------------------------------------------------
 

----------------------------
***      References      ***
----------------------------

 [1] M. Newman & M. Girvan. Finding and evaluating community structure in networks. Physical Review, E 69 (026113), 2004.

 [2] M. Newman. Finding community structure in networks using the eigenvectors of matrices. Physical Review E, 74(3):036104, 2006.
