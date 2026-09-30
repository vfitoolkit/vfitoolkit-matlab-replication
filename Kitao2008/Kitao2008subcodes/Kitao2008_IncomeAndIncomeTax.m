function [TaxI,TaxI_nonlinearpart,I]=Kitao2008_IncomeAndIncomeTax(e,a,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k,tau_E1)
% Taxable income and income tax liability of a household, for each of the three tax systems
% considered by Kitao (2008). Kept in one place so that the household budget constraint (ReturnFn,
% ConsumptionFn) and the government revenue (TaxFn, IncomeTaxFn) can never drift apart.
%
% Outputs: TaxI                is the whole income tax liability
%          TaxI_nonlinearpart  is just the tau_a0*() progressive part of that liability
%          I                   is total taxable income
%
% TaxI_nonlinearpart is reported separately because Kitao (2008), Section 3.2, calibrates the
% parameter tau_a2 so that "the share of the government expenditures raised by the non-linear part
% of the function equals 65%". So the 65% target in Table 2 is the non-linear part as a share of
% total tax revenue, and NOT the whole income tax as a share of total tax revenue.
% (Note that under taxincome=2 and taxincome=3 the flat taxes tau_k and tau_E1 are, like tau_I,
% part of the 'rest', not part of the non-linear part.)
%
% taxincome=1: tax total income I                        (benchmark)
% taxincome=2: tax capital income ra at flat rate tau_k, and the rest of income Itilde=I-ra
%              using the progressive schedule plus tau_I  (Section 5.1 of Kitao (2008))
% taxincome=3: tax entrepreneurial business income IE1 at flat rate tau_E1, and individual
%              income (IE2=ra for entrepreneurs, I for workers) using the progressive
%              schedule plus tau_I                       (Section 5.2 of Kitao (2008))
%
% The if I>0 style guards just keep the progressive part at zero on a non-positive tax base
% (at the optimum these never bind, but they avoid complex numbers if a general eqm iteration
% ever tries a negative r).

TaxI=0; % just to make GPU happy
TaxI_nonlinearpart=0;

ra=r*a; % IE2 in Kitao's notation for taxincome=3: Section 5.2 defines it as "the capital income earned on assets a"

% IK is the base of the flat capital income tax of taxincome=2, and is NOT the same as ra. Section
% 5.1 taxes "the return from riskless saving", so for an entrepreneur it covers only the part of the
% assets that is not invested in the business: "If only part of his assets are invested, i.e., k<=a,
% the remaining (a-k) earns a riskless return which is added to the tax base of the entrepreneur as
% capital income". When k>a the entrepreneur is a net borrower and so has no riskless saving at all.
% This is what makes entrepreneurial investment fall as tau_k falls, as in Fig 4(e) of Kitao (2008):
% a low tau_k makes saving attractive relative to business investment, because business profit sits
% in Itilde and is taxed at the higher tau_I. With IK=ra the entrepreneur would get the tau_k break
% on the whole asset stock regardless of where it is deployed, that margin would disappear, and
% entrepreneurial capital would move the other way.
IK=0; % just to make GPU happy

if e==0 % Workers
    I=w*eta+ra; % Eqn (4) of Kitao (2008)
    IE1=0; % not used for workers, but must be assigned to keep the GPU happy
    IK=ra; % all of a worker's assets are riskless saving
else % Entrepreneurs. Their production problem is static
    [k,n,output]=Kitao2008_StaticEntrepreneurProblem(a,eta,r,w,d,phi,delta,theta,upsilon1,upsilon2,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);

    rbar=r;
    if k>a
        rbar=r+phi; % Cost of borrowing
    end

    I=output-delta*k-rbar*(k-a)-w*(n-eta)*(n>eta); % Eqn (8) of Kitao (2008). Note: (n-eta)*(n>eta) is just max(n-eta,0)
    % Note: the -rbar*(k-a) is because
    % "If the agent is a net borrower, i.e., k>a, the interest payment for the borrowing is deducted
    % as operational costs. If only part of his assets are invested, i.e., k<=a, the remaining (a-k)
    % earns a riskless return which is added to the tax base of the entrepreneur as capital income."

    IK=r*(a-k)*(k<a); % only the assets not invested in the business earn the riskless return

    IE1=output-delta*k-r*k-phi*(k-a)*(k>a)-w*(n-eta)*(n>eta); % Eqn (11) of Kitao (2008), only used when taxincome=3
    % Note: I has -rbar*(k-a), while IE1 has -r*k-phi*(k-a)*(k>a)
    % When k>a, rbar=r+phi. When k<=a, rbar=r
end

if taxincome==1 % This is baseline model where income is taxed
    if I>0
        TaxI_nonlinearpart=tau_a0*(I-(tau_a2+I^(-tau_a1))^(-1/tau_a1));
        TaxI=TaxI_nonlinearpart+tau_I*I;
    end
elseif taxincome==2 % Seperate taxation of labor income and capital income
    Itilde=I-IK;
    TaxI=tau_k*IK;
    if Itilde>0
        TaxI_nonlinearpart=tau_a0*(Itilde-(tau_a2+Itilde^(-tau_a1))^(-1/tau_a1));
        TaxI=TaxI+TaxI_nonlinearpart+tau_I*Itilde;
    end
elseif taxincome==3 % Tax entrepreneur business income
    if e==0 % Workers are taxed on their income exactly as in the benchmark
        if I>0
            TaxI_nonlinearpart=tau_a0*(I-(tau_a2+I^(-tau_a1))^(-1/tau_a1));
            TaxI=TaxI_nonlinearpart+tau_I*I;
        end
    else % Entrepreneurs: flat tax on business income, progressive on individual (capital) income
        if IE1>0
            TaxI=tau_E1*IE1;
        end
        if ra>0
            TaxI_nonlinearpart=tau_a0*(ra-(tau_a2+ra^(-tau_a1))^(-1/tau_a1));
            TaxI=TaxI+TaxI_nonlinearpart+tau_I*ra;
        end
    end
end

end
