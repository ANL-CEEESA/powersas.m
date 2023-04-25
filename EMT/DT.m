function [t_new,x_new,y_new,LTE,LTE2,LTE3] = DT(t,x,f,GLO,h,K)    
% DT performs one-step integration of differential transformation
% method to solve differential equation dx/dt = f(t,x,GLO);
%
% Input
%   k      current order, 0,1,2,...
%   x      current states
%   f      differential equation
%   GLO    parameters in the differential equation
%   h      time step
%   K      The degree of power series: X(0)+X(1)*h+...X(K)*h^K
%
% Output
%   t_new  next time instant
%   x_new  next states of state variables
%   y_new  next states of algebraic variables
%   LTE    local truncation error
%
%  See also rk4, Parareal
% 


    %% Function handle to compute intermidate variables 
    fun_y   = str2func(strcat(func2str(f),'_AE'));
    fun_idx = str2func(strcat(func2str(f),'_idx'));

    %% Function handle for DT model X(k+1)=FB(X(0:k)) of differential equation dx/dt = f(x)   
    F = str2func(strcat(func2str(f),'_DT'));
    
    %% compute intermidate variables y
    idx     = fun_idx(GLO);
    y       = fun_y(x,idx,GLO);
    
    %% Define matrix X Y to store coefficients, Xh to store terms
    Nx  = length(x);
    X = zeros(Nx,K+1); %%why NAN
    
    Ny = length(y);
    Y = zeros(Ny,K+1); 

    Xh= zeros(Nx,K+1);
    Yh= zeros(Ny,K+1);
    %% Initialize the 0th order coefficients
    X(:,1) = x; 
    Y(:,1) = y; 

    Xh(:,1) = x;   % Term(0) = X(0)*h^0
    Yh(:,1) = y;
    %% Perform recursive calculation of high order coefficients and terms
    for k = 0:K-1 
%     for k = 0:K-2 
        %% DT Equation: X(k+1)=F(X(0:k)) or XX(k+2)=F(XX(1:k+1))
        [X,Y] = F(k,X,Y,GLO,idx); 
        Xh(:,k+2) = X(:,k+2)*h^(k+1); % Term(k+1) = X(k+1)*h^(k+1) or XX(k+2)*h^(k+1)
        Yh(:,k+2) = Y(:,k+2)*h^(k+1); % Term(k+1) = X(k+1)*h^(k+1) or XX(k+2)*h^(k+1)

    end
    
    %% sum over all columns of matrix Xh to obtain x_new as a column vector
    x_new = sum(Xh(:,1:end-1),2);
   
    %% 
%     y_new = sum(Yh,2); 
    y_new = fun_y(x_new,idx,GLO);
    %% Update time
    t_new = t + h;
    eta = 0.0000000000000000001;
    %% Local truncation error: last term; maximum value over all variables
    % khuang: we can use relative error, since different variables have
    % differnet dynamics and quantities
    infcof   = max(abs(Xh(:,end)./(x_new+eta)));
    trc   = infcof^(1/K);
    LTE   = infcof/(1-trc);
%     LTE   = infcof;
%     LTE   = max(abs(Xh(:,end)./x_new));
    LTE3   = max(abs(Xh(:,end-1)./(x_new+eta)));

    % relative error used to change order
    LTE2 = infcof/min([norm(Xh(:,end-1)./Xh(:,end),Inf), sqrt(norm(Xh(:,end-2)./Xh(:,end),Inf)), sqrt(norm(Xh(:,end-3)./Xh(:,end-1),Inf))]);
%     LTE2   = min(abs(Xh(:,end)./abs(Xh(:,end-1))))
end