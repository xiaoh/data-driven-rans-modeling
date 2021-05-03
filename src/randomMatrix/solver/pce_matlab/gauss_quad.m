function [xi,w] = gauss_quad(m)
%% Gauss quadratures for several orthogonal polynomials
%{
---------------------------------------------------------------------------
Created by:                          Date:              Comment:
Felipe Uribe                         Oct 2014           PC exercises
furibec@unal.edu.co
Universidad Nacional de Colombia
Manizales Campus
---------------------------------------------------------------------------
Comments:
input:   m     - quadrature points
         type  - chosen orthogonal polynomial weighting function:
                 'he_prob': Hermite probabilist
                 'he_phys': Hermite physicist
                 'legen'  : Legendre
                 'cheby'  : Chebyshev
                 'lague'  : Laguerre
output:  xi - orthogonal polynomial roots
         w  - weights
---------------------------------------------------------------------------
EXAMPLES:
1. For the probabilist Hermite polynomials there is no much information.
Therefore, by comparing with the moments of a standard Gaussian variable 
(http://en.wikipedia.org/wiki/Normal_distribution#Moments), one can see 
that the implementation is correct
---
clear; 
sigma = 1;
df    = @(n) prod(n:-2:1);                   % double factorial
mom   = @(p) (sigma^p)*df(p-1)*mod(p+1,2);   % p-th moment for sigma = 1
m     = 10;
[xi, w] = quad_GH(m,'he_prob');   % probabilist Hermite
for i = 1:(2*m-1)
   [i sum(xi.^i.*w)  mom(i)]
end

2. For the rest of orthogonal polynomials, compare results with tables in:
http://mathworld.wolfram.com/Hermite-GaussQuadrature.html
http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
http://mathworld.wolfram.com/Chebyshev-GaussQuadrature.html
http://mathworld.wolfram.com/Laguerre-GaussQuadrature.html
---------------------------------------------------------------------------
Based on:
1. Press et al. - "Numerical recipes"
---------------------------------------------------------------------------
%}
% syms x;
% switch type
%    case 'he_prob'   % probabilist Hermite polynomials
% probabilist Hermite polynomials

syms x;
P    = cell(m+1,1);
P{1} = 1;       % H_1 = 1
P{2} = [1 0];   % H_2 = x
for n = 2:m
    P{n+1} = [P{n} 0] - (n-1)*[0 0 P{n-1}];   % recursive formula
end
% weight function
w_x = exp(-x^2/2)/sqrt(2*pi);
% int limits
a = -Inf;   b = Inf;
      
%    case 'he_phys'   % physicist Hermite polynomials
%       P    = cell(m+1,1);
%       P{1} = 1;       % H_1 = 1
%       P{2} = [2 0];   % H_2 = 2x
%       for n = 2:m
%          P{n+1} = 2*[P{n} 0] - (2*(n-1))*[0 0 P{n-1}];   % recursive formula
%       end
%       % weight function
%       w_x = exp(-x^2);
%       % int limits
%       a = -Inf;   b = Inf;
%       
%    case 'legen'   % Legendre polynomials
%       P    = cell(m+1,1);
%       P{1} = 1;       % P_0 = 1
%       P{2} = [1 0];   % P_1 = x
%       for n = 2:m
%          P{n+1} = ((2*n-1)*[P{n} 0] - (n-1)*[0 0 P{n-1}])/n; 
%       end
%       % weight function
%       w_x = 1;
%       % int limits
%       a = -1;   b = 1;
%       
%    case 'cheby'   % Chebyshev polynomials (First kind)
%       P    = cell(m+1,1);
%       P{1} = 1;       % T_0 = 1
%       P{2} = [1 0];   % T_1 = x
%       for n = 2:m
%          P{n+1} = 2*[P{n} 0] - [0 0 P{n-1}]; 
%       end
%       % weight function
%       w_x = 1/sqrt(1-x^2);
%       % int limits
%       a = -1;   b = 1;
%       
%    case 'lague'   % generalized Laguerre polynomials (set alpha manually here)
%       P     = cell(m+1,1);
%       alpha = 0;              % alpha = 0, simple Laguerre polys
%       P{1}  = 1;              % L_0 = 1
%       P{2}  = [-1 1+alpha];   % L_1 = 1+alpha -x
%       for n = 2:m
%          P{n+1} = (((2*n-1+alpha)*[0 P{n}]-[P{n} 0]) - (n-1+alpha)*[0 0 P{n-1}])/n; 
%       end
%       % weight function
%       w_x   = (x^alpha)*exp(-x);
%       % int limits
%       a = 0;   b = Inf;
%       
%    case 'jacobi'   % Jacobi polynomials
%       alpha = param(1);
%       beta  = param(2);
%       P     = cell(m+1,1);
%       P{1}  = 1;       % P_0 = 1
%       P{2}  = (1/2)*[alpha+beta+2 alpha-beta];   
%       for n = 2:m
%          P{n+1} = ((2*n-1)*[P{n} 0] - (n-1)*[0 0 P{n-1}])/n; 
%       end
%       % weight function
%       w_x = ((1-x)^alpha)*((1+x)^beta);
%       % int limits
%       a = -1;   b = 1;
%       
%    otherwise
%       error('You must choose between "he_prob", "he_phys", "legen", "cheby", "lague"');
% end

%% find the roots
xi = sort(roots(P{m+1}));
      
%% compute the weights using general formula 
% http://en.wikipedia.org/wiki/Gaussian_quadrature#General_formula_for_the_weights
w  = zeros(m,1);
id = eye(m);
%
warning('off','all');  % Turn off warning in polyfit when using high hermite
for i = 1:m
   lag_pol = poly2sym(polyfit(xi,id(:,i),m-1),x);   % lagrange polynomial
   w(i)    = int(lag_pol*w_x, x, a, b);
end  
%
if ~isreal(w)   % ~max m = 66
   error('m is too large. The weights cannot be complex.');
end

%% plots
%{
% zeros and weights
close all;
figure;
plot(xi,w,'bo','LineWidth',2); 
set(gca,'FontSize',12);
grid minor; axis tight; ylim([min(w)-0.1 max(w)+0.1]);
xlabel('Zeros','FontSize',12);
ylabel('Weights','FontSize',12);

% polynomials
xx = min(xi):0.01:max(xi);
yy = zeros(m,length(xx));
figure;
for i = 1:m
   yy(i,:) = polyval(P{i},xx);
end
h = plot(xx,yy,'LineWidth',2);
legend(h,'Location','Best'); 
set(gca,'FontSize',12);
grid minor; axis tight; 
% ylim([-10 20]);   % for Laguerre 
xlabel('x','FontSize',12);
ylabel('P_n(x)','FontSize',12);
%}

return;
%%END