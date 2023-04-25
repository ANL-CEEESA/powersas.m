 
clc; clear; close;
addpath('EMT')

testSystem = 'twoarea';

[t_final, sol_final] = calEMT(testSystem);

plotEMTResults();

