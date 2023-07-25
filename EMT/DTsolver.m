function [t_dtm,sol_dtm]=DTsolver(TimeRange,Resolution,x,N_DTM,Sys)
Gen_num=length(Sys.GenIdx);
%% Asign space to store the series
DTA=zeros(Gen_num,N_DTM+1); % machine number, term number, dimention of each term
THETA=DTA;OMG=DTA; LAMBDA_F=DTA;LAMBDA_D=DTA;LAMBDA_Q1=DTA;LAMBDA_Q2=DTA;
V1 =DTA; EFD=DTA; P1= DTA; P2= DTA; PM=DTA; W_R=DTA;
SIN=DTA; COS=DTA; ID =DTA; IQ =DTA; 
LAMBDA_AD=DTA;  LAMBDA_AQ=DTA;  PE=DTA;  VD_PP=DTA;  VQ_PP=DTA;
VPP_A=DTA;  VPP_B=DTA;  VPP_C=DTA;  
VT=DTA;  VTP=DTA;

VBUS5 = zeros(3,N_DTM+1);
VBUS6 = VBUS5;  VBUS7 =VBUS5; VBUS8 =VBUS5; VBUS9 =VBUS5; VBUS10 =VBUS5;  VBUS11 =VBUS5; 
ILINE1 = VBUS5;  ILINE2 = VBUS5;  ILINE3 = VBUS5;  ILINE4 = VBUS5; ILINE5 = VBUS5; ILINE6 = VBUS5;
ILINE7 = VBUS5;  ILINE8 = VBUS5;  ILINE9 = VBUS5;  ILINE10= VBUS5; ILINE11= VBUS5; ILINE12= VBUS5;
ILOAD1 = VBUS5;  ILOAD2 = VBUS5;  

%% Coefficients of t^0, initial value
DTA(:,1)=x(1:4,1); THETA(:,1)=x(5:8,1); OMG(:,1)=x(9:12,1);
LAMBDA_F(:,1)=x(13:16,1); LAMBDA_D(:,1)=x(17:20,1); LAMBDA_Q1(:,1)=x(21:24,1);  LAMBDA_Q2(:,1)=x(25:28,1); V1(:,1)=x(29:32,1); 
efd  = x(33:36,1);  efd = min(max(efd,Sys.Emin),Sys.Emax);   EFD(:,1)= efd; 
p1   = x(37:40,1);  p1  = min(max(p1,Sys.GVmin),Sys.GVmax);  P1(:,1) = p1;
P2(:,1)= x(41:44,1);
PM(:,1)= P2(:,1)-Sys.Gdt.*OMG(:,1);
W_R(:,1)   =OMG(:,1)+1;
VBUS5(:,1) =x(45:47,1);  VBUS6(:,1) =x(48:50,1);   VBUS7(:,1) =x(51:53,1);   VBUS8(:,1) =x(54:56,1);  
VBUS9(:,1) =x(57:59,1);  VBUS10(:,1)=x(60:62,1);   VBUS11(:,1)=x(63:65,1);  
ILINE1(:,1)=x(66:68,1);  ILINE2(:,1)=x(69:71,1);   ILINE3(:,1)=x(72:74,1);   ILINE4(:,1)=x(75:77,1);
ILINE5(:,1)=x(78:80,1);  ILINE6(:,1)=x(81:83,1);   ILINE7(:,1)=x(84:86,1);   ILINE8(:,1)=x(87:89,1);
ILINE9(:,1)=x(90:92,1);  ILINE10(:,1)=x(93:95,1);  ILINE11(:,1)=x(96:98,1);  ILINE12(:,1)=x(99:101,1);
ILOAD1(:,1)=x(102:104,1);  ILOAD2(:,1)=x(105:107,1);
VTP(:,1)   =x(108:111,1);

% some other terms
SIN(:,1) =sin(THETA(:,1));  % sin(theta)  no need to add '
COS(:,1) =cos(THETA(:,1));  % cos(theta)

ID(:,1) = 2/3*( [ILINE1(1,1); ILINE2(1,1); ILINE3(1,1); ILINE4(1,1)].*COS(:,1)   +  [ILINE1(2,1); ILINE2(2,1); ILINE3(2,1); ILINE4(2,1)].*( COS(:,1)*cos(2*pi/3)+SIN(:,1)*sin(2*pi/3) )  +  [ILINE1(3,1); ILINE2(3,1); ILINE3(3,1); ILINE4(3,1)].*( COS(:,1)*cos(2*pi/3)-SIN(:,1)*sin(2*pi/3) )  );
IQ(:,1) = 2/3*( -[ILINE1(1,1); ILINE2(1,1); ILINE3(1,1); ILINE4(1,1)].*SIN(:,1)  -  [ILINE1(2,1); ILINE2(2,1); ILINE3(2,1); ILINE4(2,1)].*( SIN(:,1)*cos(2*pi/3)-COS(:,1)*sin(2*pi/3) )  -  [ILINE1(3,1); ILINE2(3,1); ILINE3(3,1); ILINE4(3,1)].*( SIN(:,1)*cos(2*pi/3)+COS(:,1)*sin(2*pi/3) )  );

LAMBDA_AD(:,1)  = Sys.Lad_pp.*(-ID(:,1)+LAMBDA_F(:,1)./Sys.Lfd+LAMBDA_D(:,1)./Sys.L1d);
LAMBDA_AQ(:,1)  = Sys.Laq_pp.*(-IQ(:,1)+LAMBDA_Q1(:,1)./Sys.L1q+LAMBDA_Q2(:,1)./Sys.L2q);  

PE(:,1)    = LAMBDA_AD(:,1).*IQ(:,1)-LAMBDA_AQ(:,1).*ID(:,1); 

VD_PP(:,1)=Sys.Vd_pp_coeff1.*ID(:,1)+Sys.Vd_pp_coeff2.*LAMBDA_F(:,1)+Sys.Vd_pp_coeff3.*LAMBDA_D(:,1)...
           -W_R(:,1).*Sys.Laq_pp.*( LAMBDA_Q1(:,1)./Sys.L1q + LAMBDA_Q2(:,1)./Sys.L2q ) + Sys.Lad_pp./Sys.Lfd.*EFD(:,1) ;

VQ_PP(:,1)=Sys.Vq_pp_coeff1.*IQ(:,1)+Sys.Vq_pp_coeff2.*LAMBDA_Q1(:,1)+Sys.Vq_pp_coeff3.*LAMBDA_Q2(:,1)...
           +W_R(:,1).*Sys.Lad_pp.*( LAMBDA_F(:,1)./Sys.Lfd + LAMBDA_D(:,1)./Sys.L1d );

VPP_A(:,1) = COS(:,1).*VD_PP(:,1)                                  -SIN(:,1).*VQ_PP(:,1);
VPP_B(:,1) = ( COS(:,1)*cos(2*pi/3)+SIN(:,1)*sin(2*pi/3) ).*VD_PP(:,1)  -( SIN(:,1)*cos(2*pi/3)-COS(:,1)*sin(2*pi/3) ).*VQ_PP(:,1);
VPP_C(:,1) = ( COS(:,1)*cos(2*pi/3)-SIN(:,1)*sin(2*pi/3) ).*VD_PP(:,1)  -( SIN(:,1)*cos(2*pi/3)+COS(:,1)*sin(2*pi/3) ).*VQ_PP(:,1);
VT2(:,1)  =  abs(VD_PP(:,1)+1i*VQ_PP(:,1)).^2;
VT(:,1)   =  abs(VD_PP(:,1)+1i*VQ_PP(:,1));

%% preparation 
SimuTime = TimeRange(2)-TimeRange(1);
t_dtm=(0:Resolution:SimuTime)'+TimeRange(1);
% hh = TimeRange(2)-t_dtm(end);
% if abs(hh)>1e-5
%   t_dtm = [t_dtm;TimeRange(2)];
% end

%% store the results at each time step
% assign space
res_dta   = zeros(Gen_num,size(t_dtm,1));
res_theta = res_dta;
res_omg   = res_dta; 
res_lambda_f  =  res_dta;  
res_lambda_d  =  res_dta;  
res_lambda_q1 =  res_dta;  
res_lambda_q2 =  res_dta;  
res_v1    = res_dta;   
res_efd   = res_dta;
res_p1    = res_dta;
res_p2    = res_dta;
res_vbus5  = zeros(3,size(t_dtm,1));  res_vbus6  = res_vbus5;   res_vbus7  = res_vbus5;   
res_vbus8  = res_vbus5;   res_vbus9  = res_vbus5;   res_vbus10 = res_vbus5;   res_vbus11 = res_vbus5;
res_iline1 = res_vbus5;   res_iline2 = res_vbus5;   res_iline3 = res_vbus5;   res_iline4 = res_vbus5;
res_iline5 = res_vbus5;   res_iline6 = res_vbus5;   res_iline7 = res_vbus5;   res_iline8 = res_vbus5;
res_iline9 = res_vbus5;   res_iline10= res_vbus5;   res_iline11= res_vbus5;   res_iline12= res_vbus5;
res_iload1 = res_vbus5;   res_iload2 = res_vbus5;
res_vtp    = res_dta;

% store the first one
res_dta(:,1)  =DTA(:,1);
res_theta(:,1)=THETA(:,1);
res_omg(:,1)  =OMG(:,1,1);
res_lambda_f(:,1)  =LAMBDA_F(:,1); 
res_lambda_d(:,1)  =LAMBDA_D(:,1); 
res_lambda_q1(:,1) =LAMBDA_Q1(:,1); 
res_lambda_q2(:,1) =LAMBDA_Q2(:,1); 
res_v1(:,1)  =V1(:,1); 
res_efd(:,1) =EFD(:,1);
res_p1(:,1)  =P1(:,1);
res_p2(:,1)  =P2(:,1);
res_vbus5(:,1) = VBUS5(:,1);   res_vbus6(:,1)= VBUS6(:,1);    res_vbus7(:,1)= VBUS7(:,1);   res_vbus8(:,1)= VBUS8(:,1);  
res_vbus9(:,1) = VBUS9(:,1);   res_vbus10(:,1)= VBUS10(:,1);  res_vbus11(:,1) = VBUS11(:,1);
res_iline1(:,1) = ILINE1(:,1); res_iline2(:,1) = ILINE2(:,1);   res_iline3(:,1) = ILINE3(:,1);  res_iline4(:,1) = ILINE4(:,1);
res_iline5(:,1) = ILINE5(:,1); res_iline6(:,1) = ILINE6(:,1);   res_iline7(:,1) = ILINE7(:,1);  res_iline8(:,1) = ILINE8(:,1);
res_iline9(:,1) = ILINE9(:,1); res_iline10(:,1)= ILINE10(:,1);  res_iline11(:,1)= ILINE11(:,1); res_iline12(:,1)= ILINE12(:,1);
res_iload1(:,1) = ILOAD1(:,1); res_iload2(:,1) = ILOAD2(:,1);
res_vtp(:,1)    = VTP(:,1);

%% begin the for loop to calculate 
Tvec=zeros(1,N_DTM+1);
for Numstep=1:size(t_dtm,1)-1  % loop of time step
  hhh = t_dtm(Numstep+1)-t_dtm(Numstep);
  Tvec(1:N_DTM+1)=hhh.^(0:N_DTM);
  idx = zeros(N_DTM+1,1); idx(1) =1; % for constant terms
% below are used to check limits
def_dt  = ( V1(:,1).*Sys.Ek-EFD(:,1).*Sys.Lad./Sys.Rfd )./Sys.Ete.*Sys.Rfd./Sys.Lad;
dP1_dt  = ( (Sys.Pref-OMG(:,1))./Sys.Gr-P1(:,1) )./Sys.Gt1; 
for k = 1:N_DTM     % loop of order

% synchronous generator differential equations
DTA(:,k+1)   = 1/k*2*pi*60*OMG(:,k);
THETA(:,k+1) = 1/k*2*pi*60*(OMG(:,k)+1*idx(k));
OMG(:,k+1)   = 1/k*Sys.w0*(PM(:,k) - PE(:,k) - Sys.D.*OMG(:,k)./Sys.w0)./Sys.H/2;
LAMBDA_F(:,k+1)   = 1/k*120*pi*(EFD(:,k)-Sys.Rfd./Sys.Lfd.*(LAMBDA_F(:,k)-LAMBDA_AD(:,k)));
LAMBDA_D(:,k+1)   = 1/k*(-120)*pi*Sys.R1d./Sys.L1d.*(LAMBDA_D(:,k) -LAMBDA_AD(:,k));
LAMBDA_Q1(:,k+1)  = 1/k*(-120)*pi*Sys.R1q./Sys.L1q.*(LAMBDA_Q1(:,k)-LAMBDA_AQ(:,k));
LAMBDA_Q2(:,k+1)  = 1/k*(-120)*pi*Sys.R2q./Sys.L2q.*(LAMBDA_Q2(:,k)-LAMBDA_AQ(:,k));
V1(:,k+1)         = 1/k*( -VTP(:,k)+Sys.Vref*idx(k)-V1(:,k)-Sys.Eta.*(VT(:,k)-VTP(:,k))./Sys.T_vt )./Sys.Etb;
EFD_high   =  1/k*( V1(:,k).*Sys.Ek-EFD(:,k).*Sys.Lad./Sys.Rfd )./Sys.Ete.*Sys.Rfd./Sys.Lad;
EFD_high(def_dt(EFD(:,1)>=Sys.Emax)>0)=0;
EFD_high(def_dt(EFD(:,1)<=Sys.Emin)<0)=0;
EFD(:,k+1) =  EFD_high;

P1_high = 1/k*( (Sys.Pref*idx(k)-OMG(:,k))./Sys.Gr-P1(:,k) )./Sys.Gt1; 
P1_high(dP1_dt(P1(:,1)>=Sys.GVmax)>0)=0;
P1_high(dP1_dt(P1(:,1)<=Sys.GVmin)<0)=0;
P1(:,k+1) = P1_high;

% P2(:,k+1)= 1/k*(P1(:,k)-P2(:,k)+P1(:,k+1).*Sys.Gt2)./Sys.Gt3;
P2(:,k+1)= 1/k*(P1(:,k)-P2(:,k)+k*P1(:,k+1).*Sys.Gt2)./Sys.Gt3;

% network differential equations
VBUS5(:,k+1) = 1/k*Sys.Cnode_3Phase(:,:,5)*(ILINE1(:,k)-ILINE5(:,k));
VBUS6(:,k+1) = 1/k*Sys.Cnode_3Phase(:,:,6)*(ILINE2(:,k)+ILINE5(:,k)-ILINE6(:,k));
VBUS7(:,k+1) = 1/k*Sys.Cnode_3Phase(:,:,7)*(ILINE6(:,k)-ILINE7(:,k)-ILINE8(:,k)-ILOAD1(:,k));
VBUS8(:,k+1) = 1/k*Sys.Cnode_3Phase(:,:,8)*(ILINE7(:,k)+ILINE8(:,k)-ILINE9(:,k)-ILINE10(:,k));
VBUS9(:,k+1) = 1/k*Sys.Cnode_3Phase(:,:,9)*(ILINE9(:,k)+ILINE10(:,k)-ILINE11(:,k)-ILOAD2(:,k));
VBUS10(:,k+1)= 1/k*Sys.Cnode_3Phase(:,:,10)*(ILINE4(:,k)+ILINE11(:,k)-ILINE12(:,k));
VBUS11(:,k+1)= 1/k*Sys.Cnode_3Phase(:,:,11)*(ILINE3(:,k)+ILINE12(:,k));

ILINE1(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,1)*( [VPP_A(1,k); VPP_B(1,k); VPP_C(1,k)]- VBUS5(:,k)  -Sys.Gen_R(:,:,1)*ILINE1(:,k) );  % IBR voltage
ILINE2(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,2)*( [VPP_A(2,k); VPP_B(2,k); VPP_C(2,k)]- VBUS6(:,k)  -Sys.Gen_R(:,:,2)*ILINE2(:,k) ); 
ILINE3(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,3)*( [VPP_A(3,k); VPP_B(3,k); VPP_C(3,k)]- VBUS11(:,k) -Sys.Gen_R(:,:,3)*ILINE3(:,k) ); 
ILINE4(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,4)*( [VPP_A(4,k); VPP_B(4,k); VPP_C(4,k)]- VBUS10(:,k) -Sys.Gen_R(:,:,4)*ILINE4(:,k) ); 

ILINE5(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,5)*(  VBUS5(:,k)-VBUS6(:,k)-Sys.Rline_3Phase(:,:,5) * ILINE5(:,k)  );
ILINE6(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,6)*(  VBUS6(:,k)-VBUS7(:,k)-Sys.Rline_3Phase(:,:,6) * ILINE6(:,k)  );
ILINE7(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,7)*(  VBUS7(:,k)-VBUS8(:,k)-Sys.Rline_3Phase(:,:,7) * ILINE7(:,k)  );
ILINE8(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,8)*(  VBUS7(:,k)-VBUS8(:,k)-Sys.Rline_3Phase(:,:,8) * ILINE8(:,k)  );
ILINE9(:,k+1)= 1/k*Sys.Lline_3Phase(:,:,9)*(  VBUS8(:,k)-VBUS9(:,k)-Sys.Rline_3Phase(:,:,9) * ILINE9(:,k)  );
ILINE10(:,k+1)=1/k*Sys.Lline_3Phase(:,:,10)*( VBUS8(:,k)-VBUS9(:,k)-Sys.Rline_3Phase(:,:,10) * ILINE10(:,k)  );
ILINE11(:,k+1)=1/k*Sys.Lline_3Phase(:,:,11)*( VBUS9(:,k)-VBUS10(:,k)-Sys.Rline_3Phase(:,:,11) * ILINE11(:,k)  );
ILINE12(:,k+1)=1/k*Sys.Lline_3Phase(:,:,12)*( VBUS10(:,k)-VBUS11(:,k)-Sys.Rline_3Phase(:,:,12) * ILINE12(:,k)  );

ILOAD1(:,k+1)=1/k*Sys.Lload_3Phase(:,:,1)*(  VBUS7(:,k)-Sys.Rload_3Phase(:,:,1) * ILOAD1(:,k)  );
ILOAD2(:,k+1)=1/k*Sys.Lload_3Phase(:,:,2)*(  VBUS9(:,k)-Sys.Rload_3Phase(:,:,2) * ILOAD2(:,k)  );

VTP(:,k+1)=1/k*(VT(:,k)-VTP(:,k))./Sys.T_vt;

% below is algebric part
PM(:,k+1)= P2(:,k+1)-Sys.Gdt.*OMG(:,k+1);
W_R(:,k+1)   = OMG(:,k+1)+1*idx(k+1);

COS(:,k+1)= -1/(k)*sum(repmat([k:-1:1],Gen_num,1).*SIN(:,1:k).*THETA(:,k+1:-1:2),2);
SIN(:,k+1)=  1/(k)*sum(repmat([k:-1:1],Gen_num,1).*COS(:,1:k).*THETA(:,k+1:-1:2),2);

ID(:,k+1) = 2/3*sum(( [ILINE1(1,1:k+1); ILINE2(1,1:k+1); ILINE3(1,1:k+1); ILINE4(1,1:k+1)].*COS(:,k+1:-1:1)  +  [ILINE1(2,1:k+1); ILINE2(2,1:k+1); ILINE3(2,1:k+1); ILINE4(2,1:k+1)].*( COS(:,k+1:-1:1)*cos(2*pi/3)+SIN(:,k+1:-1:1)*sin(2*pi/3) )  +  [ILINE1(3,1:k+1); ILINE2(3,1:k+1); ILINE3(3,1:k+1); ILINE4(3,1:k+1)].*( COS(:,k+1:-1:1)*cos(2*pi/3)-SIN(:,k+1:-1:1)*sin(2*pi/3) )  ),2);
IQ(:,k+1) = 2/3*sum((-[ILINE1(1,1:k+1); ILINE2(1,1:k+1); ILINE3(1,1:k+1); ILINE4(1,1:k+1)].*SIN(:,k+1:-1:1)  -  [ILINE1(2,1:k+1); ILINE2(2,1:k+1); ILINE3(2,1:k+1); ILINE4(2,1:k+1)].*( SIN(:,k+1:-1:1)*cos(2*pi/3)-COS(:,k+1:-1:1)*sin(2*pi/3) )  -  [ILINE1(3,1:k+1); ILINE2(3,1:k+1); ILINE3(3,1:k+1); ILINE4(3,1:k+1)].*( SIN(:,k+1:-1:1)*cos(2*pi/3)+COS(:,k+1:-1:1)*sin(2*pi/3) )  ),2);

LAMBDA_AD(:,k+1)  = Sys.Lad_pp.*(-ID(:,k+1)+LAMBDA_F(:,k+1)./Sys.Lfd+LAMBDA_D(:,k+1)./Sys.L1d);
LAMBDA_AQ(:,k+1)  = Sys.Laq_pp.*(-IQ(:,k+1)+LAMBDA_Q1(:,k+1)./Sys.L1q+LAMBDA_Q2(:,k+1)./Sys.L2q);  

PE(:,k+1)    = sum(LAMBDA_AD(:,1:k+1).*IQ(:,k+1:-1:1)-LAMBDA_AQ(:,1:k+1).*ID(:,k+1:-1:1),2); 

VD_PP(:,k+1)=Sys.Vd_pp_coeff1.*ID(:,k+1)+Sys.Vd_pp_coeff2.*LAMBDA_F(:,k+1)+Sys.Vd_pp_coeff3.*LAMBDA_D(:,k+1)...
           -sum(W_R(:,1:k+1).*Sys.Laq_pp.*( LAMBDA_Q1(:,k+1:-1:1)./Sys.L1q + LAMBDA_Q2(:,k+1:-1:1)./Sys.L2q ),2) + Sys.Lad_pp./Sys.Lfd.*EFD(:,k+1) ;

VQ_PP(:,k+1)=Sys.Vq_pp_coeff1.*IQ(:,k+1)+Sys.Vq_pp_coeff2.*LAMBDA_Q1(:,k+1)+Sys.Vq_pp_coeff3.*LAMBDA_Q2(:,k+1)...
           +sum(W_R(:,1:k+1).*Sys.Lad_pp.*( LAMBDA_F(:,k+1:-1:1)./Sys.Lfd + LAMBDA_D(:,k+1:-1:1)./Sys.L1d ),2);

VPP_A(:,k+1) = sum(COS(:,1:k+1).*VD_PP(:,k+1:-1:1)                                  -SIN(:,1:k+1).*VQ_PP(:,k+1:-1:1),2);
VPP_B(:,k+1) = sum(( COS(:,1:k+1)*cos(2*pi/3)+SIN(:,1:k+1)*sin(2*pi/3) ).*VD_PP(:,k+1:-1:1)  -( SIN(:,1:k+1)*cos(2*pi/3)-COS(:,1:k+1)*sin(2*pi/3) ).*VQ_PP(:,k+1:-1:1),2);
VPP_C(:,k+1) = sum(( COS(:,1:k+1)*cos(2*pi/3)-SIN(:,1:k+1)*sin(2*pi/3) ).*VD_PP(:,k+1:-1:1)  -( SIN(:,1:k+1)*cos(2*pi/3)+COS(:,1:k+1)*sin(2*pi/3) ).*VQ_PP(:,k+1:-1:1),2);

VT2(:,k+1) = sum(VD_PP(:,1:k+1).*VD_PP(:,k+1:-1:1)+ VQ_PP(:,1:k+1).*VQ_PP(:,k+1:-1:1),2); 
VT(:,k+1)  = (VT2(:,k+1) - sum(VT(:,2:k).*VT(:,k:-1:2),2))./(2*VT(:,1));

end

%%  Evaluation, contain only differential equations
dta_2=sum(DTA.*repmat(Tvec,[Gen_num,1]),2);
theta_2=sum(THETA.*repmat(Tvec,[Gen_num,1]),2);
omg_2=sum(OMG.*repmat(Tvec,[Gen_num,1]),2);
lambda_f_2=sum(LAMBDA_F.*repmat(Tvec,[Gen_num,1]),2);
lambda_d_2=sum(LAMBDA_D.*repmat(Tvec,[Gen_num,1]),2);
lambda_q1_2=sum(LAMBDA_Q1.*repmat(Tvec,[Gen_num,1]),2);
lambda_q2_2=sum(LAMBDA_Q2.*repmat(Tvec,[Gen_num,1]),2); 
vtp_2=sum(VTP.*repmat(Tvec,[Gen_num,1]),2); 
v1_2=sum(V1.*repmat(Tvec,[Gen_num,1]),2); 
efd_2=sum(EFD.*repmat(Tvec,[Gen_num,1]),2); 
p1_2=sum(P1.*repmat(Tvec,[Gen_num,1]),2); 
p2_2=sum(P2.*repmat(Tvec,[Gen_num,1]),2); 
vbus5_2=sum(VBUS5.*repmat(Tvec,[3,1]),2);
vbus6_2=sum(VBUS6.*repmat(Tvec,[3,1]),2);
vbus7_2=sum(VBUS7.*repmat(Tvec,[3,1]),2);
vbus8_2=sum(VBUS8.*repmat(Tvec,[3,1]),2);
vbus9_2=sum(VBUS9.*repmat(Tvec,[3,1]),2);
vbus10_2=sum(VBUS10.*repmat(Tvec,[3,1]),2);
vbus11_2=sum(VBUS11.*repmat(Tvec,[3,1]),2);
iline1_2=sum(ILINE1.*repmat(Tvec,[3,1]),2);
iline2_2=sum(ILINE2.*repmat(Tvec,[3,1]),2);
iline3_2=sum(ILINE3.*repmat(Tvec,[3,1]),2);
iline4_2=sum(ILINE4.*repmat(Tvec,[3,1]),2);
iline5_2=sum(ILINE5.*repmat(Tvec,[3,1]),2);
iline6_2=sum(ILINE6.*repmat(Tvec,[3,1]),2);
iline7_2=sum(ILINE7.*repmat(Tvec,[3,1]),2);
iline8_2=sum(ILINE8.*repmat(Tvec,[3,1]),2);
iline9_2=sum(ILINE9.*repmat(Tvec,[3,1]),2);
iline10_2=sum(ILINE10.*repmat(Tvec,[3,1]),2);
iline11_2=sum(ILINE11.*repmat(Tvec,[3,1]),2);
iline12_2=sum(ILINE12.*repmat(Tvec,[3,1]),2);
iload1_2=sum(ILOAD1.*repmat(Tvec,[3,1]),2);
iload2_2=sum(ILOAD2.*repmat(Tvec,[3,1]),2);
 
% Save the results    
res_dta(:,Numstep+1)=dta_2;
res_theta(:,Numstep+1)=theta_2;
res_omg(:,Numstep+1)=omg_2;
res_lambda_f(:,Numstep+1)=lambda_f_2;
res_lambda_d(:,Numstep+1)=lambda_d_2;
res_lambda_q1(:,Numstep+1)=lambda_q1_2;
res_lambda_q2(:,Numstep+1)=lambda_q2_2;
res_v1(:,Numstep+1)=v1_2;
res_efd(:,Numstep+1)=efd_2;
res_p1(:,Numstep+1)=p1_2;
res_p2(:,Numstep+1)=p2_2;
res_vtp(:,Numstep+1)=vtp_2;
res_vbus5(:,Numstep+1)=vbus5_2;    
res_vbus6(:,Numstep+1)=vbus6_2; 
res_vbus7(:,Numstep+1)=vbus7_2;    
res_vbus8(:,Numstep+1)=vbus8_2;
res_vbus9(:,Numstep+1)=vbus9_2;    
res_vbus10(:,Numstep+1)=vbus10_2;
res_vbus11(:,Numstep+1)=vbus11_2;    
res_iline1(:,Numstep+1)=iline1_2; 
res_iline2(:,Numstep+1)=iline2_2;
res_iline3(:,Numstep+1)=iline3_2; 
res_iline4(:,Numstep+1)=iline4_2;
res_iline5(:,Numstep+1)=iline5_2; 
res_iline6(:,Numstep+1)=iline6_2;
res_iline7(:,Numstep+1)=iline7_2; 
res_iline8(:,Numstep+1)=iline8_2;
res_iline9(:,Numstep+1)=iline9_2; 
res_iline10(:,Numstep+1)=iline10_2;
res_iline11(:,Numstep+1)=iline11_2; 
res_iline12(:,Numstep+1)=iline12_2;
res_iload1(:,Numstep+1)=iload1_2;
res_iload2(:,Numstep+1)=iload2_2;

%%  Update Initial value
% differential part
DTA(:,1)=dta_2;
THETA(:,1)=theta_2;
OMG(:,1)=omg_2;
LAMBDA_F(:,1)=lambda_f_2;
LAMBDA_D(:,1)=lambda_d_2;
LAMBDA_Q1(:,1)=lambda_q1_2;
LAMBDA_Q2(:,1)=lambda_q2_2;
V1(:,1)=v1_2;
EFD(:,1)=efd_2;
P1(:,1)=p1_2;
P2(:,1)=p2_2;
VTP(:,1)=vtp_2;
VBUS5(:,1)=vbus5_2;
VBUS6(:,1)=vbus6_2;
VBUS7(:,1)=vbus7_2;
VBUS8(:,1)=vbus8_2;
VBUS9(:,1)=vbus9_2;
VBUS10(:,1)=vbus10_2;
VBUS11(:,1)=vbus11_2;
ILINE1(:,1)=iline1_2;
ILINE2(:,1)=iline2_2;
ILINE3(:,1)=iline3_2;
ILINE4(:,1)=iline4_2;
ILINE5(:,1)=iline5_2;
ILINE6(:,1)=iline6_2;
ILINE7(:,1)=iline7_2;
ILINE8(:,1)=iline8_2;
ILINE9(:,1)=iline9_2;
ILINE10(:,1)=iline10_2;
ILINE11(:,1)=iline11_2;
ILINE12(:,1)=iline12_2;
ILOAD1(:,1)=iload1_2;
ILOAD2(:,1)=iload2_2;

% algebric part 
PM(:,1)= P2(:,1)-Sys.Gdt.*OMG(:,1);
W_R(:,1) = OMG(:,1)+1;
SIN(:,1) =sin(THETA(:,1));  % sin(theta)  no need to add '
COS(:,1) =cos(THETA(:,1));  % cos(theta)

ID(:,1) = 2/3*( [ILINE1(1,1); ILINE2(1,1); ILINE3(1,1); ILINE4(1,1)].*COS(:,1)   +  [ILINE1(2,1); ILINE2(2,1); ILINE3(2,1); ILINE4(2,1)].*( COS(:,1)*cos(2*pi/3)+SIN(:,1)*sin(2*pi/3) )  +  [ILINE1(3,1); ILINE2(3,1); ILINE3(3,1); ILINE4(3,1)].*( COS(:,1)*cos(2*pi/3)-SIN(:,1)*sin(2*pi/3) )  );
IQ(:,1) = 2/3*( -[ILINE1(1,1); ILINE2(1,1); ILINE3(1,1); ILINE4(1,1)].*SIN(:,1)  -  [ILINE1(2,1); ILINE2(2,1); ILINE3(2,1); ILINE4(2,1)].*( SIN(:,1)*cos(2*pi/3)-COS(:,1)*sin(2*pi/3) )  -  [ILINE1(3,1); ILINE2(3,1); ILINE3(3,1); ILINE4(3,1)].*( SIN(:,1)*cos(2*pi/3)+COS(:,1)*sin(2*pi/3) )  );

LAMBDA_AD(:,1)  = Sys.Lad_pp.*(-ID(:,1)+LAMBDA_F(:,1)./Sys.Lfd+LAMBDA_D(:,1)./Sys.L1d);
LAMBDA_AQ(:,1)  = Sys.Laq_pp.*(-IQ(:,1)+LAMBDA_Q1(:,1)./Sys.L1q+LAMBDA_Q2(:,1)./Sys.L2q);  

PE(:,1)    = LAMBDA_AD(:,1).*IQ(:,1)-LAMBDA_AQ(:,1).*ID(:,1); 

VD_PP(:,1)=Sys.Vd_pp_coeff1.*ID(:,1)+Sys.Vd_pp_coeff2.*LAMBDA_F(:,1)+Sys.Vd_pp_coeff3.*LAMBDA_D(:,1)...
           -W_R(:,1).*Sys.Laq_pp.*( LAMBDA_Q1(:,1)./Sys.L1q + LAMBDA_Q2(:,1)./Sys.L2q ) + Sys.Lad_pp./Sys.Lfd.*EFD(:,1) ;

VQ_PP(:,1)=Sys.Vq_pp_coeff1.*IQ(:,1)+Sys.Vq_pp_coeff2.*LAMBDA_Q1(:,1)+Sys.Vq_pp_coeff3.*LAMBDA_Q2(:,1)...
           +W_R(:,1).*Sys.Lad_pp.*( LAMBDA_F(:,1)./Sys.Lfd + LAMBDA_D(:,1)./Sys.L1d );

VPP_A(:,1) = COS(:,1).*VD_PP(:,1)                                  -SIN(:,1).*VQ_PP(:,1);
VPP_B(:,1) = ( COS(:,1)*cos(2*pi/3)+SIN(:,1)*sin(2*pi/3) ).*VD_PP(:,1)  -( SIN(:,1)*cos(2*pi/3)-COS(:,1)*sin(2*pi/3) ).*VQ_PP(:,1);
VPP_C(:,1) = ( COS(:,1)*cos(2*pi/3)-SIN(:,1)*sin(2*pi/3) ).*VD_PP(:,1)  -( SIN(:,1)*cos(2*pi/3)+COS(:,1)*sin(2*pi/3) ).*VQ_PP(:,1);

VT2(:,1)  =  abs(VD_PP(:,1)+1i*VQ_PP(:,1)).^2;
VT(:,1)   =  abs(VD_PP(:,1)+1i*VQ_PP(:,1));

end

sol_dtm = zeros(110,size(t_dtm,1));  % generator state variables 
sol_dtm(1:4,:) = res_dta;
sol_dtm(5:8,:) = res_theta;
sol_dtm(9:12,:) = res_omg;
sol_dtm(13:16,:) = res_lambda_f;
sol_dtm(17:20,:) = res_lambda_d;
sol_dtm(21:24,:) = res_lambda_q1;
sol_dtm(25:28,:) = res_lambda_q2;
sol_dtm(29:32,:) = res_v1;
sol_dtm(33:36,:) = res_efd;
sol_dtm(37:40,:) = res_p1;
sol_dtm(41:44,:) = res_p2;

sol_dtm(45:47,:) = res_vbus5;
sol_dtm(48:50,:) = res_vbus6;
sol_dtm(51:53,:) = res_vbus7;
sol_dtm(54:56,:) = res_vbus8;
sol_dtm(57:59,:) = res_vbus9;
sol_dtm(60:62,:) = res_vbus10;
sol_dtm(63:65,:) = res_vbus11;
sol_dtm(66:68,:) = res_iline1;
sol_dtm(69:71,:) = res_iline2;
sol_dtm(72:74,:) = res_iline3;
sol_dtm(75:77,:) = res_iline4;
sol_dtm(78:80,:) = res_iline5;
sol_dtm(81:83,:) = res_iline6;
sol_dtm(84:86,:) = res_iline7;
sol_dtm(87:89,:) = res_iline8;
sol_dtm(90:92,:) = res_iline9;
sol_dtm(93:95,:) = res_iline10; 
sol_dtm(96:98,:) = res_iline11;
sol_dtm(99:101,:) = res_iline12;
sol_dtm(102:104,:) = res_iload1;
sol_dtm(105:107,:) = res_iload2;
sol_dtm(108:111,:) = res_vtp;

