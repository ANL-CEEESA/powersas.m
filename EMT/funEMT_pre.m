function [xp] = funEMT_pre(t,x,Sys)
xp    = zeros(size(x));
% global GLO
dta   = x(1:4,1);
theta = x(4+1:2*4,1);
omg   = x(2*4+1:3*4,1);
Lambda_F   = x(3*4+1:4*4,1);
Lambda_D   = x(4*4+1:5*4,1);
Lambda_Q1  = x(5*4+1:6*4,1);
Lambda_Q2  = x(6*4+1:7*4,1);

V1    = x(7*4+1:8*4,1);
efd   = x(8*4+1:9*4,1);
% efd   = min(max(efd,Sys.Emin),Sys.Emax);

P1    = x(9*4+1:10*4,1);
% P1    = min(max(P1,Sys.GVmin),Sys.GVmax);

P2    = x(10*4+1:11*4,1);
Pm    = P2-Sys.Gdt.*omg;
Vbus  = x(11*4+1:11*4+3*(11-4),1);
Iline = x(11*4+3*(11-4)+1:11*4+3*(11-4)+3*12,1);
Iload = x(11*4+3*(11-4)+3*12+1:11*4+3*(11-4)+3*12+3*2,1);
Vtp   = x(107+1:107+4,1);

w_r   = omg+1;
COS =  cos(theta);
SIN =  sin(theta);

Id = 2/3*( Iline(1:3:3*4).*COS  + Iline(2:3:3*4).*( COS*cos(2*pi/3)+SIN*sin(2*pi/3) ) + Iline(3:3:3*4).*( COS*cos(2*pi/3)-SIN*sin(2*pi/3) )  );
Iq = 2/3*( -Iline(1:3:3*4).*SIN - Iline(2:3:3*4).*( SIN*cos(2*pi/3)-COS*sin(2*pi/3) ) - Iline(3:3:3*4).*( SIN*cos(2*pi/3)+COS*sin(2*pi/3) )  );

Lambda_ad  = Sys.Lad_pp.*(-Id+Lambda_F./Sys.Lfd+Lambda_D./Sys.L1d);
Lambda_aq  = Sys.Laq_pp.*(-Iq+Lambda_Q1./Sys.L1q+Lambda_Q2./Sys.L2q);  
Pe    = (Lambda_ad.*Iq-Lambda_aq.*Id); 

Vd_pp=Sys.Vd_pp_coeff1.*Id+Sys.Vd_pp_coeff2.*Lambda_F+Sys.Vd_pp_coeff3.*Lambda_D...
      -w_r.*Sys.Laq_pp.*( Lambda_Q1./Sys.L1q + Lambda_Q2./Sys.L2q ) + Sys.Lad_pp./Sys.Lfd.*efd ;

Vq_pp=Sys.Vq_pp_coeff1.*Iq+Sys.Vq_pp_coeff2.*Lambda_Q1+Sys.Vq_pp_coeff3.*Lambda_Q2...
      +w_r.*Sys.Lad_pp.*( Lambda_F./Sys.Lfd + Lambda_D./Sys.L1d );

%  Inv_Park(:,:,k)*[0; Vd_pp(k); Vq_pp(k)]
Vpp_a   = COS.*Vd_pp                                  -SIN.*Vq_pp;
Vpp_b   = ( COS*cos(2*pi/3)+SIN*sin(2*pi/3) ).*Vd_pp  -( SIN*cos(2*pi/3)-COS*sin(2*pi/3) ).*Vq_pp;
Vpp_c   = ( COS*cos(2*pi/3)-SIN*sin(2*pi/3) ).*Vd_pp  -( SIN*cos(2*pi/3)+COS*sin(2*pi/3) ).*Vq_pp;

Vt   =  sqrt(Vd_pp.^2+Vq_pp.^2);  % in fact it is internal voltage

% ODEs
xp(1:4,1)  = 2*pi*60*omg;
xp(4+1 :  2*4,1)  = 2*pi*60*(omg+1);
xp(2*4+1: 3*4,1)  = Sys.w0*(Pm - Pe - Sys.D.*omg./Sys.w0)./Sys.H/2;
xp(3*4+1: 4*4,1)  = 120*pi*(efd-Sys.Rfd./Sys.Lfd.*(Lambda_F-Lambda_ad));
xp(4*4+1: 5*4,1)  = -120*pi*Sys.R1d./Sys.L1d.*(Lambda_D-Lambda_ad);
xp(5*4+1: 6*4,1)  = -120*pi*Sys.R1q./Sys.L1q.*(Lambda_Q1-Lambda_aq);
xp(6*4+1: 7*4,1)  = -120*pi*Sys.R2q./Sys.L2q.*(Lambda_Q2-Lambda_aq);

xp(7*4+1: 8*4,1)  =  ( -Vtp+Sys.Vref-V1-Sys.Eta.*(Vt-Vtp)./Sys.T_vt )./Sys.Etb; % substitute d(Vt')=(Vt-Vtp)./Sys.T_vt
xp(8*4+1: 9*4,1)  =  ( V1.*Sys.Ek-efd.*Sys.Lad./Sys.Rfd )./Sys.Ete.*Sys.Rfd./Sys.Lad;  
 
xp(9*4+1: 10*4,1)  = ( (Sys.Pref-omg)./Sys.Gr-P1 )./Sys.Gt1; 
xp(10*4+1: 11*4,1) = (P1-P2+( (Sys.Pref-omg)./Sys.Gr-P1 )./Sys.Gt1.*Sys.Gt2)./Sys.Gt3;

xp(44+1:44+3,1) = Sys.Cnode_3Phase(:,:,5)*([sum(Iline(([1 ]-1)*3+1,1));sum(Iline(([1 ]-1)*3+2,1));sum(Iline(([1 ]-1)*3+3,1))]-[sum(Iline(([5 ]-1)*3+1,1));sum(Iline(([5 ]-1)*3+2,1));sum(Iline(([5 ]-1)*3+3,1))]);
xp(47+1:47+3,1) = Sys.Cnode_3Phase(:,:,6)*([sum(Iline(([2 5 ]-1)*3+1,1));sum(Iline(([2 5 ]-1)*3+2,1));sum(Iline(([2 5 ]-1)*3+3,1))]-[sum(Iline(([6 ]-1)*3+1,1));sum(Iline(([6 ]-1)*3+2,1));sum(Iline(([6 ]-1)*3+3,1))]);
xp(50+1:50+3,1) = Sys.Cnode_3Phase(:,:,7)*([sum(Iline(([6 ]-1)*3+1,1));sum(Iline(([6 ]-1)*3+2,1));sum(Iline(([6 ]-1)*3+3,1))]-[sum(Iline(([7 8 ]-1)*3+1,1));sum(Iline(([7 8 ]-1)*3+2,1));sum(Iline(([7 8 ]-1)*3+3,1))]-[sum(Iload(([1 ]-1)*3+1,1));sum(Iload(([1 ]-1)*3+2,1));sum(Iload(([1 ]-1)*3+3,1))]);
xp(53+1:53+3,1) = Sys.Cnode_3Phase(:,:,8)*([sum(Iline(([7 8 ]-1)*3+1,1));sum(Iline(([7 8 ]-1)*3+2,1));sum(Iline(([7 8 ]-1)*3+3,1))]-[sum(Iline(([9 10 ]-1)*3+1,1));sum(Iline(([9 10 ]-1)*3+2,1));sum(Iline(([9 10 ]-1)*3+3,1))]);
xp(56+1:56+3,1) = Sys.Cnode_3Phase(:,:,9)*([sum(Iline(([9 10 ]-1)*3+1,1));sum(Iline(([9 10 ]-1)*3+2,1));sum(Iline(([9 10 ]-1)*3+3,1))]-[sum(Iline(([11 ]-1)*3+1,1));sum(Iline(([11 ]-1)*3+2,1));sum(Iline(([11 ]-1)*3+3,1))]-[sum(Iload(([2 ]-1)*3+1,1));sum(Iload(([2 ]-1)*3+2,1));sum(Iload(([2 ]-1)*3+3,1))]);
xp(59+1:59+3,1) = Sys.Cnode_3Phase(:,:,10)*([sum(Iline(([4 11 ]-1)*3+1,1));sum(Iline(([4 11 ]-1)*3+2,1));sum(Iline(([4 11 ]-1)*3+3,1))]-[sum(Iline(([12 ]-1)*3+1,1));sum(Iline(([12 ]-1)*3+2,1));sum(Iline(([12 ]-1)*3+3,1))]);
xp(62+1:62+3,1) = Sys.Cnode_3Phase(:,:,11)*([sum(Iline(([3 12 ]-1)*3+1,1));sum(Iline(([3 12 ]-1)*3+2,1));sum(Iline(([3 12 ]-1)*3+3,1))]);

xp(65+1:65+3,1)= Sys.Lline_3Phase(:,:,1)*( [Vpp_a(1); Vpp_b(1); Vpp_c(1)]- Vbus((5-5)*3+1:(5-5)*3+3)  -Sys.Gen_R(:,:,1)*Iline((1-1)*3+1:(1-1)*3+3));  % IBR voltage
xp(68+1:68+3,1)= Sys.Lline_3Phase(:,:,2)*( [Vpp_a(2); Vpp_b(2); Vpp_c(2)]- Vbus((6-5)*3+1:(6-5)*3+3)  -Sys.Gen_R(:,:,2)*Iline((2-1)*3+1:(2-1)*3+3)); 
xp(71+1:71+3,1)= Sys.Lline_3Phase(:,:,3)*( [Vpp_a(3); Vpp_b(3); Vpp_c(3)]- Vbus((11-5)*3+1:(11-5)*3+3)-Sys.Gen_R(:,:,3)*Iline((3-1)*3+1:(3-1)*3+3)); 
xp(74+1:74+3,1)= Sys.Lline_3Phase(:,:,4)*( [Vpp_a(4); Vpp_b(4); Vpp_c(4)]- Vbus((10-5)*3+1:(10-5)*3+3)-Sys.Gen_R(:,:,4)*Iline((4-1)*3+1:(4-1)*3+3)); 

xp(77+1:77+3,1)=Sys.Lline_3Phase(:,:,5)*(  Vbus((5-5)*3+1:(5-5)*3+3)-Vbus((6-5)*3+1:(6-5)*3+3)-Sys.Rline_3Phase(:,:,5) * Iline((5-1)*3+1:(5-1)*3+3)  );
xp(80+1:80+3,1)=Sys.Lline_3Phase(:,:,6)*(  Vbus((6-5)*3+1:(6-5)*3+3)-Vbus((7-5)*3+1:(7-5)*3+3)-Sys.Rline_3Phase(:,:,6) * Iline((6-1)*3+1:(6-1)*3+3)  );
xp(83+1:83+3,1)=Sys.Lline_3Phase(:,:,7)*(  Vbus((7-5)*3+1:(7-5)*3+3)-Vbus((8-5)*3+1:(8-5)*3+3)-Sys.Rline_3Phase(:,:,7) * Iline((7-1)*3+1:(7-1)*3+3)  );
xp(86+1:86+3,1)=Sys.Lline_3Phase(:,:,8)*(  Vbus((7-5)*3+1:(7-5)*3+3)-Vbus((8-5)*3+1:(8-5)*3+3)-Sys.Rline_3Phase(:,:,8) * Iline((8-1)*3+1:(8-1)*3+3)  );
xp(89+1:89+3,1)=Sys.Lline_3Phase(:,:,9)*(  Vbus((8-5)*3+1:(8-5)*3+3)-Vbus((9-5)*3+1:(9-5)*3+3)-Sys.Rline_3Phase(:,:,9) * Iline((9-1)*3+1:(9-1)*3+3)  );
xp(92+1:92+3,1)=Sys.Lline_3Phase(:,:,10)*(  Vbus((8-5)*3+1:(8-5)*3+3)-Vbus((9-5)*3+1:(9-5)*3+3)-Sys.Rline_3Phase(:,:,10) * Iline((10-1)*3+1:(10-1)*3+3)  );
xp(95+1:95+3,1)=Sys.Lline_3Phase(:,:,11)*(  Vbus((9-5)*3+1:(9-5)*3+3)-Vbus((10-5)*3+1:(10-5)*3+3)-Sys.Rline_3Phase(:,:,11) * Iline((11-1)*3+1:(11-1)*3+3)  );
xp(98+1:98+3,1)=Sys.Lline_3Phase(:,:,12)*(  Vbus((10-5)*3+1:(10-5)*3+3)-Vbus((11-5)*3+1:(11-5)*3+3)-Sys.Rline_3Phase(:,:,12) * Iline((12-1)*3+1:(12-1)*3+3)  );

xp(101+1:101+3,1)=Sys.Lload_3Phase(:,:,1)*(  Vbus((7-5)*3+1:(7-5)*3+3)-Sys.Rload_3Phase(:,:,1) * Iload((1-1)*3+1:(1-1)*3+3)  );
xp(104+1:104+3,1)=Sys.Lload_3Phase(:,:,2)*(  Vbus((9-5)*3+1:(9-5)*3+3)-Sys.Rload_3Phase(:,:,2) * Iload((2-1)*3+1:(2-1)*3+3)  );

xp(107+1:107+4,1)=(Vt-Vtp)./Sys.T_vt;

% GLO.V1 =[GLO.V1,  Vt(1)];
% GLO.V2 =[GLO.V2,  Vt(2)];
% GLO.V3 =[GLO.V3,  Vt(3)];
% GLO.V4 =[GLO.V4,  Vt(4)];
% plot(0:2*10^(-6):t,GLO.V2(1:4:end))
end
