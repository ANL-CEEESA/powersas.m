clc; clear; close;
addpath('EMT')

res = runPowerSAS('emt','twoarea');

plotEMTResults(res);

