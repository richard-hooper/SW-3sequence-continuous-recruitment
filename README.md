This repository contains supplementary files to the following article:

Efficient designs for three-sequence stepped wedge trials with continuous recruitment

Richard Hooper(1), Olivier Quintin(1), Jessica Kasza(2)
(1) Queen Mary University of London, London, UK (2) Monash University, Melbourne, Australia

Address for correspondence:
Richard Hooper,
Queen Mary University of London,
Wolfson Institute of Population Health, 
Yvonne Carter Building,
58 Turner Street,
London E1 2AB

r.l.hooper@qmul.ac.uk

[This article is awaiting publication - citation details will follow]


The repository comprises:

1. Mata and Stata code (see the file sw_3seq_cr.do) for calculating and plotting the variance
   of the treatment effect estimator from a three-sequence stepped wedge trial with continuous
   recruitment. Mata is a C-style language which runs in Stata (Stata Corp, College Station,
   Texas USA).

3. R code for a simulation study covering the scenarios considered in the PATHWEIGH example.
   This code generates empirical power for the scenarios in Table 1 of the paper. There are
   three files: ThreestepSWcode.R contains the R functions needed for the simulation;
   ThreestepSWscenarios.R contains the code needed to run the simulation study; and
   scenarios_power.csv is a spreadsheet containing the parameters for the considered scenarios,
   along with the theoretical power associated with each.
