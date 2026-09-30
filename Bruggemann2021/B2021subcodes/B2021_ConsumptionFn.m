function c=B2021_ConsumptionFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,lumpsum,pension,tau_c,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)
% B2021_ConsumptionFn Household consumption implied by current choices.
% Returns consumption for workers, entrepreneurs, and retirees. For young
% entrepreneurs, firm-side objects are recovered from B2021_EntrepreneurStaticProblem using
% the fixed own-labor input lbar.

c=0;

if e==0 && age==1 % Young workers' problem
    I=w*l*eta+r*a; % income, y^w

    % Calculate the income tax
    TaxI=B2021_IncomeTaxFn(I,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);
        
    c=(I-TaxI+lumpsum+a-aprime)/(1+tau_c);
    
elseif e==1 && age==1 % Entrepreneurs' problem

    [k,n,output]=B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar);
    
    I=output-delta*k-r*(k-a)-w*n; % Income, y^e [note: -w*n is equal to -w*nbar+w*lbar
        
    % Calculate the income tax
    TaxI=B2021_IncomeTaxFn(I,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);
    
    c=(I-TaxI+lumpsum+a-aprime)/(1+tau_c);
        
elseif age==2 % Retiree's problem (for both e=0 & e=1)
    I=r*a+pension; % income, y^r

    % Calculate the income tax
    TaxI=B2021_IncomeTaxFn(I,d,tau_s,ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6);
        
    c=(I-TaxI+lumpsum+a-aprime)/(1+tau_c);
    
end %end function












end
