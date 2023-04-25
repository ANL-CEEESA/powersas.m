function [t,x,y,hh,KK,bx,by,bz] = multiwindow(daeModelEquations,simu_T,x0,GLO,method,solverOptions,bx0)
format long 
% multiwindow performs general multi-stage numerical integration. It
% applies to RK4, DT or other use-specified numerical solvers
% 
% Inputs
%       daeModelEquations  the struct for equations related to a DAE model
%       simu_T             the simulation window: start, step length, end
%       x0                 initial values of the DAE
%       GLO                parameters in the DAE model
%       method             solution method
%       solverOptions      settings and parameters of the solution method
% Outputs
%       t                  vector of time steps in the whole simulation
%       x                  trajectory of state variables
%       y                  trajectory of algebraic variables
%       hh                 the time step used in each simulation step
%       KK                 the order used in each simulation step
% 
% See also Parareal, rk4, DT, getDAEModelEquations, getDAESolver

f = daeModelEquations.fun_DAE;
fun_idx = daeModelEquations.fun_idx;

%%
K = solverOptions.K;
% FS1_VS0 = solverOptions.FS1_VS0; 
% FO1_VO0 = solverOptions.FO1_VO0;
hmax = solverOptions.hmax;
hmin = solverOptions.hmin;
Kmax = solverOptions.Kmax;
Kmin = solverOptions.Kmin;
% eps1 = solverOptions.eps1;
% eps2 = solverOptions.eps2;
% eps3 = solverOptions.eps3;
% eps4 = solverOptions.eps4;
tol =  solverOptions.RBtol;
% q1 = solverOptions.q1;
% q2 = solverOptions.q2;
% dK = solverOptions.dK;
%%
a = simu_T(1); % start time of simulation
b = simu_T(2); % end time of simulation
h = simu_T(3); % time step length

%% Judge if x0 is column vector or not. If rox vector, convert to column vector.
m = size(x0,1);  % number of rows
if m == 1
   x0 = x0';
   m  =size(x0,1);
end
bm = size(bx0,1);  % number of rows
if bm == 1
   bx0 = bx0';
end
%%
fun_y   = str2func(strcat(func2str(f),'_AE'));
idx     = fun_idx(GLO);
y0 = fun_y(x0,idx,GLO);
by0 = fun_y(bx0,idx,GLO);
%% Initial values of t and x
% mite = ceil((b-a)./hmin);
mite = 20000;
bt = zeros(1,mite);
bx = zeros(m,mite);
by = zeros(size(y0,1),mite);
t  = zeros(1,mite);
x  = zeros(m,mite);
y  = zeros(size(y0,1),mite);
hh = zeros(1,mite);
KK = zeros(1,mite);
t(1) = a;
x(:,1) = x0;     %initial conditions
bt(:,1) = a;
bx(:,1) = bx0;

y(:,1) = y0;
by(:,1) = by0;
hh(:,1)= h;
KK(:,1)= K;

%% Implement the integration method
ratio_b4= nan;
theta_max=2;



ite = 1 ;

% K    = -(1/2)*log(tol);
M = 1; 
fac1= 0.98;
fac2 = 1;

% fac1 = 0;
% fac2 = 0;
% 

p= 1 ;
pf1 = 2*GLO.machine;
pf2 = 2*GLO.machine+4*GLO.machine*GLO.machine+idx.Nx*(idx.Nx+idx.Ny)+(idx.Nx);
pf3 = (idx.Nx);
% pf1 = 17*GLO.machine/2;
% pf2 = 27*GLO.machine/2+4*GLO.machine*GLO.machine+idx.Nx*(idx.Nx+idx.Ny)+(idx.Nx);
% pf2 = 0;
% pf3 = 0;
while t(ite)+h <= b

switch solverOptions.solver
    case 'DT'
         
        [t_new,x_new,y_new, LTE,LTE2,LTE3] = method(t(ite),x(:,ite),f,GLO,h,K);
      ratio = LTE/h;
      K_I = 0.3/K;
      K_p = 0.4/K;
      if isnan(ratio_b4)
      else
          theta = ((tol/ratio)^(K_I))*((ratio_b4/ratio)^(K_p));

         if mod(ite,M)==0
%              Lin  = ((K+p)/(K))
%              Lin  = ((K+p)^2/(K)^2)
             Lin  = (pf1*(K+p)^2+pf2*(K+p)+pf3)/(pf1*(K)^2+pf2*(K)+pf3);
%              Lde  = ((K-p)/(K))
%              Lde  = ((K-p)^2/(K)^2)
             Lde  = (pf1*(K-p)^2+pf2*(K-p)+pf3)/(pf1*(K)^2+pf2*(K)+pf3);
             
%              trc_in =(LTE2); 
%              trc_de= (LTE3);
             
             trc_in =(LTE2)/(1-LTE2^((K+p))); 
             trc_de= (LTE3)/(1-LTE3^((K-p)));
             
             tol_in = tol;
             tol_de = tol;
             
             theta_in =1*(h*tol_in/trc_in)^(1/(K+p));
             theta_de =1*(h*tol_de/trc_de)^(1/(K-p));
             
%              
%              theta_in =1*(tol_in/trc_in)^(1/(K+p));
%              theta_de =1*(tol_de/trc_de)^(1/(K-p));

% %              
% %              theta_in =0.9*(h*tol/trc_in)^(1/(K+p));
% %              theta_de =0.9*(h*tol/trc_de)^(1/(K-p));

         if  theta> theta_max
            h_can     = theta_max*h;
         else
            h_can     = theta*h; 
         end    


                    if h_can > hmax
                        h_can = hmax;

                    elseif h < hmin
                        h_can = hmin;

                    end 

                    if  theta_in> theta_max
                       theta_in = theta_max;
                    end
                    if  theta_de> theta_max
                       theta_de = theta_max;
                    end



                    if h*theta_in > hmax
                        h_inp = hmax;

                    elseif h*theta_in < hmin
                        h_inp = hmin;
                    else
                        
                        h_inp = h*theta_in;
                    end 
                    
                    
                    if h*theta_de > hmax
                        h_dep = hmax;

                    elseif h*theta_de < hmin
                        h_dep = hmin;
                    else
                        h_dep = h*theta_de;
                    end 








                       
%              if min(hh(end-1:-1:end-M))>= h*theta&&Lin <fac1*(h_inp./h)
%              
%                         K=K+p;
%                         theta =theta_in
%              elseif Lde<fac2*(h_dep./h) &&max(hh(end-1:-1:end-M))<= h*theta
%                         K=K-p;
%                         theta =theta_de
% 
%              end

%              if (min(hh(end:-1:end-M+1))>= h_can||KK(end-1)<K)&&Lin <fac1*(h_inp./h_can)
%              
%                         K=K+p;
%                         theta =theta_in
%              elseif Lde<fac2*(h_dep./h_can) &&(max(hh(end:-1:end-M+1))<= h_can||KK(end-1)>K)
%                         K=K-p;
%                         theta =theta_de
% 
%              end


              if (Lde<fac2*(h_dep./h_can) &&(h<= h_can||(h==hmin&&KK(1,ite+1)>=K)))&&((h>= h_can||(h==hmax&&KK(1,ite+1)<=K))&&Lin <fac1*(h_inp./h_can))
                % check which one is better
                cri1 = h_inp/(pf1*(K+p)^2+pf2*(K+p)+pf3);
                cri2 = h_can/(pf1*K^2+pf2*K+pf3);
                cri3 = h_dep/(pf1*(K-p)^2+pf2*(K-p)+pf3);
                od_stratg=[cri1 cri2 cri3;
                            K+p K     K-p
                            theta_in theta theta_de
                            tol_in    tol   tol_de];
                        [~,opi] = max(od_stratg(1,:));
                        K = od_stratg(2,opi);
                        theta = od_stratg(3,opi);
                        tol = od_stratg(4,opi);
              elseif Lde<fac2*(h_dep./h_can) &&(h<= h_can||(h_can==hmin&&KK(1,ite+1)>=K))
                        K=K-p;
                        theta =theta_de;
                        tol = tol_de;

              elseif (h>= h_can||(h_can==hmax&&KK(1,ite+1)<=K))&&Lin <fac1*(h_inp./h_can)
             
                        K=K+p;
                        theta =theta_in;
                        tol = tol_in;

              end





                       
%              if Lin <fac1*(h_inp./h)
%              
%                         K=K+p;
%                         theta =theta_in
%                         tol = tol_in;
%              elseif Lde<fac2*(h_dep./h)
%                         K=K-p;
%                         theta =theta_de
%                         tol = tol_de;
% 
%              end


         end

         if  theta> theta_max
            h     = theta_max*h;
         else
            h     = theta*h; 
         end   
             
             
      end  
            ratio_b4 = ratio;        

        

    case 'RK4'
          [t_new,x_new,y_new] = method(t(ite),x(:,ite),f,GLO,h,K);

end
    
    
    
%      bhx = bx(:,ite);
%      for bht =t(ite):0.00031: t_new(end)
%          if t_new(end)-bht>=0.00031
%             [bht_new,bhx_new,bhy_new] = rk4(bht,bhx,f,GLO,0.00031,8);
%             bhx  = bhx_new(:,end);
%          else
%             [bht_new,bhx_new,bhy_new] = rk4(bht,bhx,f,GLO,t_new(end)-bht,8);
%             break
% %             bht  = t_new(end);
%          end
%      end

 if ite+1 > mite
     mite = mite + 20000;
%      bt(:,mite) = 0;
%      bx(:,mite) = 0;
%      by(:,mite) = 0;
     t(:,mite)  = 0;
     x(:,mite)  = 0;
     y(:,mite)  = 0;
     hh(:,mite) = 0;
     KK(:,mite) = 0;
 end
 
 
 
                    if h > hmax
                        h = hmax;

                    elseif h < hmin
                        h = hmin;

                    end  
             
        if K > Kmax; K = Kmax;    end
        if K < Kmin; K = Kmin;    end  
%     bt(:,ite+1) = bht_new(end);
%     bx(:,ite+1) = bhx_new(:,end);
%     by(:,ite+1) = bhy_new(:,end);
    t(:,ite+1)  = t_new;
    x(:,ite+1)  = x_new;
    y(:,ite+1)  = y_new;
    hh(:,ite+1) = h;
    KK(:,ite+1) = K; 

%       r = tol^(1/(1+K));

% theta =(tol*h/r_norm)^(1/K);
      
 
%     K;
      ite = ite+1;
    
%     if FS1_VS0 == 0
%         if LTE < eps1; h = h*q1;   end
%         if LTE > eps2; h = h*q2; end
%         if h > hmax; h = hmax;    end
%         if h < hmin; h = hmin;    end
%     end
%     
%     if FO1_VO0 == 0
%         if LTE < eps3
%             K = K-dK;   
%         end
%         if LTE > eps4 
%             K = K+dK;   
%         end
%         if K > Kmax; K = Kmax;    end
%         if K < Kmin; K = Kmin;    end    
%     end 
    if b-t(ite)<h&&b-t(ite)>0
        h = b-t(ite);
    end

end

%% transpose: use .' to avoid conjugate transpose in case some variables are complex
t  = t(:,1:ite).'; 
x  = x(:,1:ite).'; 
y  = y(:,1:ite).';
hh = hh(:,1:ite).'; 
KK = KK(:,1:ite).';
% bx = bx(:,1:ite).';
% by = by(:,1:ite).';
% bz =max(abs(x-bx),[],2);
bx = x;
by = y;
bz =max(abs(x-bx),[],2);
% bz =max((x-bx)./(max(bx+10^(-19),[],1)),[],2);
end