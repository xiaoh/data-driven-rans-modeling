%% PC applications 1D
%{
--------------------------------------------------------------------------
Created by:                       Date:           Comment:
Felipe Uribe                      Oct/2014        ---
furibec@unal.edu.co                   
Universidad Nacional de Colombia 
Manizales Campus
--------------------------------------------------------------------------
Based on:
1."Stochastic finite elements A spectral approach"
   R. Ghanem and P.D. Spanos. Rev edition 2012. Dover publications Inc.
2."Stochastic finite element methods and reliability"
   B. Sudret and A. Der Kiureghian. State of the art report.
3."Numerical methods for stochastic computations"
   D. Xiu. (2010). Princeton university press.
--------------------------------------------------------------------------
%}
function re = PCE_field(path, gammaMeanName)
    format long g;
    % Read parameters for PCE
    load([path, '/para_gamma-', gammaMeanName, '.mat']);
    %% initial parameters
    %p_order = 6;                    % order of the polynomial chaos: 0,1,2...
    xi_cdf  = @(x) normcdf(x,0,1);  % gaussian measure: homogeneous chaos
    eps_integrate = 1.0 - 1E-2;


    %% obtain the gPC Hermite expressions
    [~,PsiPol,~,PsiSqNorm,P] = Hermite_PC(p_order);  % 1D Hermite polynomials

    %% open storage for PCE coefficients and polynomials
    u_k = zeros(nCell, P);
    % convert symbolic polynomial to rows (numerical)
    num_polys = zeros(P, P);


    %% compute basis polynomials
    % Highest order of the k^{th} polynomial
    highest_order = P - [1:P] + 1;
    for k = 1:P
        num_polys(k, highest_order(k):end) = sym2poly(PsiPol{k});
    end

    %% find the PC deterministic coefficients using proyection approach
    N_quad_qt          = 6;
    [xi_q,w_q] = gauss_quad(N_quad_qt);   % Gauss-Hermite quadrature

    for cellI = 1 : nCell

        A = mean_gamma(cellI);
        B = 1;
        u_icdf = @(p) gaminv(p,A,B);   % inv cdf
        aa = 0;   
        bb = u_icdf(eps_integrate);

        for k = 1:P
            for j = 1:N_quad_qt
                u_k(cellI, k) = u_k(cellI, k) + (u_icdf(xi_cdf(xi_q(j)))*...
                    polyval(sym2poly(PsiPol{k}),xi_q(j))*w_q(j));   % Ref 3. Eq.(5.16)
            end
        end
        u_k(cellI, :) = u_k(cellI, :)./PsiSqNorm';
    end

    save([path, '/poly-', gammaMeanName, '.mat'], 'num_polys',  'u_k', 'P', 'gammaMeanName');

return;
end
% u_pdf  = @(x) gampdf(x,A,B);   % target pdf
% u_cdf  = @(x) gamcdf(x,A,B);   % cdf
% 
% %% approximate the random variable u by PCE
% N_sample       = 1e4;          % number of samples
% xi      = randn(N_sample,1);   % homogeneous chaos
% u_tilde = zeros(N_sample,1);
% for k = 1:P
%    % u_tilde = u_tilde + u_k(k)*polyval(sym2poly(PsiPol{k}),xi);
%    u_tilde = u_tilde + u_k(k)*polyval(num_polys(k, highest_order(k):end),xi);
% end

% [n1,x1]     = hist(u_tilde,ceil(sqrt(N_sample)));           % histogram
% [PC_pdf,x2] = ksdensity(u_tilde);                    % approx pdf
% [PC_cdf,x3] = ksdensity(u_tilde,'function','cdf');   % approx cdf
% 
% %% plots
% xp = aa:0.01:bb;
% set(0,'defaultTextInterpreter','latex'); 
% 
% % PDF
% figure(1); hold on;
% bar(x1,n1/trapz(x1,n1),'c','EdgeColor',[0,0,0]);
% plot(xp,u_pdf(xp),'b-','LineWidth',1);
% plot(x2,PC_pdf,'r--','LineWidth',1);  % 'k-+' 'g-.' 'r--'
% grid minor; axis tight; set(gca,'FontSize',12);
% xlabel('$$\xi$$','FontSize',15);
% ylabel('$$f_{\alpha}(\xi)$$','FontSize',15);
% legend('Hist','Exact',sprintf('%d-order PC approx',p_order),'Location','Best');
% 
% % CDF
% figure(2); 
% % cdf tail beg
% subplot(3,4,1:2); hold on;
% plot(xp,u_cdf(xp),'b-','LineWidth',2);
% plot(x3,PC_cdf,'r--','LineWidth',2); 
% grid minor; ylim([0 0.05]); set(gca,'FontSize',22);
% xlabel('$$\xi$$','FontSize',23);
% ylabel('$$F_{\alpha}(\xi)$$','FontSize',20);
% title('Tail at beginning','FontSize',25);
% % cdf tail end
% subplot(3,4,3:4); hold on;
% plot(xp,u_cdf(xp),'b-','LineWidth',2); 
% plot(x3,PC_cdf,'r--','LineWidth',2);
% grid minor; ylim([0.95 1]); set(gca,'FontSize',22);
% xlabel('$$\xi$$','FontSize',23);
% ylabel('$$F_{\alpha}(\xi)$$','FontSize',24);
% title('Tail at end','FontSize',25);
% % cdf complete
% subplot(3,4,5:12);
% plot(xp,u_cdf(xp),'b-','LineWidth',2); hold on;
% plot(x3,PC_cdf,'r--','LineWidth',2);
% grid minor; axis tight; set(gca,'FontSize',22);
% xlabel('$$\xi$$','FontSize',27);
% ylabel('$$F_{\alpha}(\xi)$$','FontSize',27);
% legend('Exact',sprintf('%d-order PC approx',p_order),'Location','Best');

%%END
