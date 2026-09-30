function F=Bruggemann2021_ReturnFn(l,e,aprime,a,age,eta,theta,r,w,sigma1,sigma2,xi,lbar,lambda,delta,gamma,upsilon,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)
% Bruggemann2021_ReturnFn computes current reward, denoted as F in toolkit notation.
%
% ACTIONS:
%   l       : Labor supply (workers only)
%   e       : Occupational choice (0 = worker, 1 = entrepreneur)
%   aprime  : Next-period assets
%
% STATES:
%   a       : Assets today
%   eta     : Worker ability
%   theta   : Entrepreneurial ability
%   age     : Age indicator (1 = young, 2 = old)
%
% These are the `always required` variables. Every input after these is a parameter.

F = -Inf;

if age==1 % Young

    if e==0 % Young worker's problem

        I = w*l*eta+r*a; % income, y^w

        % Calculate the income tax
        TaxI = B2021_IncomeTaxFn(I,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);

        % Labor supply
        labsup = l;

    else % Young entrepreneurs' problem

        % Entrepreneurs must supply the fixed labor level lbar. Keep invalid
        % l choices infeasible so the policy function continues to encode this.
        if abs(l-lbar)>1e-5
            return
        end
        
        [k,n,output] = B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar);

        I = output-delta*k-r*(k-a)-w*n; % Income, y^e

        % Calculate the income tax
        TaxI = B2021_IncomeTaxFn(I,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);

        % Labor supply: entrepreneurs must supply l = lbar
        labsup = lbar;
    end %end if e

else % Retiree's problem (for both e=0 & e=1)
    I = r*a+pension; % income, y^r

    % Calculate the income tax
    TaxI = B2021_IncomeTaxFn(I,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);

    % Labor supply: retirees do not supply labor
    labsup = 0;

end %end if age

c = (I - TaxI + lumpsum + a - aprime)/(1 + tau_c);

if c>0
    % separable utility fn
    F=(c^(1-sigma1))/(1-sigma1) - xi*(labsup^(1+sigma2))/(1+sigma2); 
end

end %end function
