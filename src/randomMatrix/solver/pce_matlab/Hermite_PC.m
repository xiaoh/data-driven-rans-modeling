function [alpha,Psi_s,Psi_p,PsiSqNorm,P] = Hermite_PC(p_order)
%% M-dimensional Hermite polynomials (homogeneous chaos)
%{
--------------------------------------------------------------------------
Created by:                       Date:           Comment:
Felipe Uribe                      Oct/2014        ---
furibec@unal.edu.co                   
Universidad Nacional de Colombia 
Manizales Campus
--------------------------------------------------------------------------
*Requires -multi_index- function
*Input:
 M           % number of K-L terms (number of random variables)
 p_order     % order of PC
--------------------------------------------------------------------------
*Output:
 alpha       % multi-index sequence
 Psi         % Hermite PC
 PsiSqNorm   % squared-norm of the Hermite PC
 P           % basis size of the Hermite PC
--------------------------------------------------------------------------
Based on:
1."Numerical methods for stochastic computations. A spectral method approach"
   D. Xiu. 2010. Princeton University Press.
2."Stochastic finite element methods and reliability"
   B. Sudret and A. Der Kiureghian. State of the art report.
--------------------------------------------------------------------------
%}

%% Calculate the basis size of Psi
P = 1; 
M = 1; % One variable only
for s = 1:p_order
   P = P + (1/factorial(s))*prod(M+(0:s-1));   % Eq.(3.13) Ref.(2)
end

%% Calculate 1D Hermite polynomials: Recurrence relation
% symbolic
syms xi;
He_s    = cell(p_order,1);
He_s{1} = sym(1);
He_s{2} = xi;
for j = 2:p_order+1
   He_s{j+1} = expand(xi*He_s{j} - (j-1)*He_s{j-1});
end
% polynomial
He_p    = cell(p_order,1);
He_p{1} = 1;       % H_1 = 1
He_p{2} = [1 0];   % H_2 = x
for n = 2:p_order+1
   He_p{n+1} = [He_p{n} 0] - (n-1)*[0 0 He_p{n-1}];   % recursive formula
end

%% Define the number of RVs
x   = cell(1,1);
H_s = cell(p_order,1);   % Hermite polinomial for 1D syms
H_p = cell(p_order,1);   % Hermite polinomial for 1D pol
for j = 1:M
   x{j} = sym(sprintf('xi_%d',j));
   for i = 1:p_order+1
      H_s{i,j} = subs(He_s{i},xi,x{j});
      H_p{i,j} = He_p{i};
   end
end

%% M-dimensional PC computation
Psi_s  = cell(P,1);   % symbolic version
Psi_p  = cell(P,1);   % polynomial version
alpha  = [0 : p_order]'; % multi_index(M,p_order);  % create the multi-index
for i = 2:P+1
    % Disable multi-dimsnsion (XH)
    % mult_s = 1;
    % mult_p = 1;
    %    for j = 1:M
    %       mult_s = mult_s*H_s{alpha(i-1,j)+1,j};
    %       mult_p = conv(mult_p,H_p{alpha(i-1,j)+1,j});
    %    end
    %     mult_s = mult_s*H_s{alpha(i-1,1)+1,1};
    %     mult_p = conv(mult_p,H_p{alpha(i-1,1)+1,1});
    %     Psi_s{i-1} = mult_s;
    %     Psi_p{i-1} = mult_p;
    Psi_s{i-1} = H_s{alpha(i-1,1)+1,1};
    Psi_p{i-1} = H_p{alpha(i-1,1)+1,1};
end

%% Calculate the square norm
PsiSqNorm  = prod(factorial(alpha),2);

%% show results
fprintf('Number of random variables (K-L terms): M = %d \n',M);
fprintf('Order of the polynomial chaos: p = %d \n',p_order);
fprintf('Total number of polynomial chaos: p = %d \n\n',P);
for k = 1:P
   fprintf('j = %d \t',k-1);
   fprintf('E[Psi^2] = %d \t\t',PsiSqNorm(k));
   fprintf('Psi   = %s \n',char(Psi_s{k}));
end

return;