function [outputArg1,outputArg2,outputArg3] = VSOO(LTE,LTE2,LTE3,SimData,h,ratio_b4)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
      tol = 10^-4;
      theta_max = 2;
      Kmax  = 30;
      Kmin  = 10;
      p     = 1;
      fac1  = 1;
      fac2  = 1;
      hmax  = SimData.maxAlpha;
      hmin  = 1*SimData.dAlpha;
      

      
      
      K = SimData.nlvl;
      ratio = LTE/h;
      K_I = 0.03/K;
      K_p = 0.04/K;
      pf1 = 1;
      pf2 = 0;
      pf3 = 0;
      if isnan(ratio_b4)
      else
             theta = ((tol/ratio)^(K_I))*((ratio_b4/ratio)^(K_p));
             Lin  = (pf1*(K+p)^2+pf2*(K+p)+pf3)/(pf1*(K)^2+pf2*(K)+pf3);
             Lde  = (pf1*(K-p)^2+pf2*(K-p)+pf3)/(pf1*(K)^2+pf2*(K)+pf3);
             
             
             trc_in =(LTE2)/(1-LTE2^((K+p))); 
             trc_de= (LTE3)/(1-LTE3^((K-p)));
             
             tol_in = tol;
             tol_de = tol;
             
             theta_in =0.9*(h*tol_in/trc_in)^(1/(K+p));
             theta_de =0.9*(h*tol_de/trc_de)^(1/(K-p));
             

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
              elseif Lde<fac2*(h_dep./h_can) &&(h<= h_can||(h_can==hmin))
                        K=K-p;
                        theta =theta_de;
                        tol = tol_de;

              elseif (h>= h_can||(h_can==hmax))&&Lin <fac1*(h_inp./h_can)
             
                        K=K+p;
                        theta =theta_in;
                        tol = tol_in;

              end
%                 cri1 = h_inp/(pf1*(K+p)^2+pf2*(K+p)+pf3);
%                 cri2 = h_can/(pf1*K^2+pf2*K+pf3);
%                 cri3 = h_dep/(pf1*(K-p)^2+pf2*(K-p)+pf3);
%                 od_stratg=[cri1 cri2 cri3;
%                             K+p K     K-p
%                             theta_in theta theta_de
%                             tol_in    tol   tol_de];
%                         [~,opi] = max(od_stratg(1,:));
%                         K = od_stratg(2,opi);
%                         theta = od_stratg(3,opi);
%                         tol = od_stratg(4,opi);
%    

         if  theta> theta_max
            h     = theta_max*h;
         else
            h     = theta*h; 
         end   
   
             
      end  
      
      if h > hmax; h = hmax;    end
      if h < hmin; h = hmin;    end 
      if K > Kmax; K = Kmax;    end
      if K < Kmin; K = Kmin;    end  
%             ratio_b4 = ratio;

outputArg1 = K;
outputArg2 = h;
outputArg3 = ratio;
end

