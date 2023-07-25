function [Sys, Laa0_pp, Lab0_pp, Laa2_pp, t_base] = fun3(Sys)

%% old, good version
t_base=1/(120*pi);

Sys.Lad=Sys.Xd-Sys.Xl;
Sys.Lfd=1./(1./(Sys.Xdp-Sys.Xl)-1./Sys.Lad);
Sys.Rfd=(Sys.Lad+Sys.Lfd)./(Sys.Td0p./t_base);
Sys.L1d=1./(1./(Sys.Xdpp-Sys.Xl)-1./Sys.Lfd-1./Sys.Lad);
Sys.R1d=(1./(1./Sys.Lad+1./Sys.Lfd)+Sys.L1d)./(Sys.Td0pp./t_base);

Sys.Laq=Sys.Xq-Sys.Xl;
Sys.L1q=1./(1./(Sys.Xqp-Sys.Xl)-1./Sys.Laq);
Sys.R1q=(Sys.Laq+Sys.L1q)./(Sys.Tq0p./t_base);  
Sys.L2q=1./(1./(Sys.Xqpp-Sys.Xl)-1./Sys.L1q-1./Sys.Laq);
Sys.R2q=(1./(1./Sys.Laq+1./Sys.L1q)+Sys.L2q)./(Sys.Tq0pp./t_base);

Sys.Lad_pp=Sys.Xdpp-Sys.Xl;
Sys.Laq_pp=Sys.Xdpp-Sys.Xl;

Laa0_pp= Sys.Xl+( Sys.Lad_pp+ Sys.Laq_pp+Sys.L0-Sys.Xl )/3;
Lab0_pp= ( 2*Sys.L0-2*Sys.Xl-Sys.Lad_pp-Sys.Laq_pp )/6;
Laa2_pp= ( Sys.Lad_pp-Sys.Laq_pp )/3;

Sys.Vd_pp_coeff1=-(Sys.Rfd./(Sys.Lfd).^2+Sys.R1d./(Sys.L1d).^2).*(Sys.Lad_pp).^2;
Sys.Vd_pp_coeff2=-(Sys.Lad_pp.*Sys.Rfd./(Sys.Lfd).^2.*(1-Sys.Lad_pp./Sys.Lfd) - Sys.Lad_pp.^2.*Sys.R1d./(Sys.L1d).^2./Sys.Lfd);
Sys.Vd_pp_coeff3=(Sys.Lad_pp.^2.*Sys.Rfd./(Sys.Lfd).^2./Sys.L1d - Sys.Lad_pp.*Sys.R1d./(Sys.L1d).^2.*(1-Sys.Lad_pp./Sys.L1d));

Sys.Vq_pp_coeff1=-(Sys.R1q./(Sys.L1q).^2+Sys.R2q./(Sys.L2q).^2).*(Sys.Laq_pp).^2;
Sys.Vq_pp_coeff2=-(Sys.Laq_pp.*Sys.R1q./(Sys.L1q).^2.*(1-Sys.Laq_pp./Sys.L1q) - Sys.Laq_pp.^2.*Sys.R2q./(Sys.L2q).^2./Sys.L1q);
Sys.Vq_pp_coeff3=(Sys.Laq_pp.^2.*Sys.R1q./(Sys.L1q).^2./Sys.L2q - Sys.Laq_pp.*Sys.R2q./(Sys.L2q).^2.*(1-Sys.Laq_pp./Sys.L2q));


