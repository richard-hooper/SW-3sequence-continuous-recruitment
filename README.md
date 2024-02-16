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

1. Mata and Stata code for calculating and plotting the variance of the treatment effect
   estimator from a three-sequence stepped wedge trial with continuous recruitment. Mata is a
   C-style language which runs in Stata (Stata Corp, College Station, Texas USA). There is one
   file (author - Richard Hooper):
   
   sw_3seq_cr.do

2. R code for the simulation study covering the scenarios considered in the PATHWEIGH example.
   This code generates empirical power for the scenarios in Table 1 of the paper. There are
   four files (author - Jessica Kasza):
   
   Scenarios.csv -- contains the parameters for each scenario;
   
   ThreeStepSWcode.R -- functions to implement the calculation of power;
   
   ThreeStepSWscenarios.R -- calculation of theoretical and empirical power;
   
   ThreeStepSWloopplot.R -- generates the nested loop plot in the paper supplement.
