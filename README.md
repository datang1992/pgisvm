pgisvm
======

This is the PGISVM code (draft version). 

This repository contains the code for the algorithm PGiSVM described in Chapter 2. There are three branches in this repository. The master branch and the pgisvm1 branch contains the code for the homogeneous classification case (i.e. Decision function: y = sign(w^T x)) while the pgisvm2 branch contains the code for the heterogeneous classication case (i.e. Decision function: y = sign(w^T x + b)). 

To run this algorithm, just use the makefile tool to compile it and run ./main. Please install Open-MP and GSL before compiling this project.

In order to change the settings for the algorithm (including the paramters and the data), please modify the source code.
