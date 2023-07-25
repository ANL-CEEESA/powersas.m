function [Lambda_F, Lambda_D, Lambda_Q1, Lambda_Q2] = calculateLambda(Sys, Id, Iq, efd)

Ifd=efd./Sys.Rfd;
Lambda_F  = -Sys.Lad.*Id+(Sys.Lfd+Sys.Lad).*Ifd;  % steady state Id1=0
Lambda_D  = -Sys.Lad.*Id+Sys.Lad.*Ifd; 
Lambda_Q1 = -Sys.Laq.*Iq; 
Lambda_Q2 = -Sys.Laq.*Iq;