function F=DFJ2010EM_ReturnFn_dgrid(d,aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)
% The endogenous medical-expense model with medical goods m as a DECISION VARIABLE on a grid (d=m), for comparison with
% DFJ2010EM_ReturnFn, which finds m from the intratemporal first-order condition by bisection. Same shifter, coinsurance,
% taxes and utility floor; here c=x-a'-q*m directly. At the floor (x<xMin) the transfer puts the person at uMin with a'=0,
% whatever d (the d is then meaningless, as the floor bundle is not on the grid).
F=-Inf;

if h==0 || h==1 % Alive (good or bad health)
    % medical-needs shifter and coinsurance rate for this health status
    mu_coeff=-Inf; s_coeff=-Inf; q=1; % placeholders (just to make gpu happy)
    if h==0 % good health
        mu_coeff=mu_coeff_healthgood; s_coeff=s_coeff_healthgood; q=q_healthgood;
    elseif h==1 % bad health
        mu_coeff=mu_coeff_healthbad; s_coeff=s_coeff_healthbad; q=q_healthbad;
    end
    mu=exp(mu_coeff+s_coeff*((zeta+xi)/normalizepsi)); % DFJ2010 eqn (21): mu(t,h,zeta,xi), C GetMu()
    IESfrac=nu/omega;
    allocfrac=q*(mu/q)^(1/omega); % intratemporal FOC c^(-nu)=(mu/q) m^(-omega): out-of-pocket q*m = allocfrac*c^IESfrac (C line 2557)
    % income and taxes (as in DeNardiFrenchJones2010_AfterTaxIncome)
    taxableincome=r*a+earnings;
    tax=0;
    if taxableincome<taxbracket1
        tax=tax+margtaxrate0*(taxableincome-0);
    else
        tax=tax+margtaxrate0*(taxbracket1-0);
    end
    if taxableincome<taxbracket2
        tax=tax+margtaxrate1*max(taxableincome-taxbracket1,0);
    else
        tax=tax+margtaxrate1*(taxbracket2-taxbracket1);
    end
    if taxableincome<taxbracket3
        tax=tax+margtaxrate2*max(taxableincome-taxbracket2,0);
    else
        tax=tax+margtaxrate2*(taxbracket3-taxbracket2);
    end
    if taxableincome<taxbracket4
        tax=tax+margtaxrate3*max(taxableincome-taxbracket3,0);
    else
        tax=tax+margtaxrate3*(taxbracket4-taxbracket3);
    end
    if taxableincome<taxbracket5
        tax=tax+margtaxrate4*max(taxableincome-taxbracket4,0);
    else
        tax=tax+margtaxrate4*(taxbracket5-taxbracket4);
    end
    if taxableincome<taxbracket6
        tax=tax+margtaxrate5*max(taxableincome-taxbracket5,0);
    else
        tax=tax+margtaxrate5*(taxbracket6-taxbracket5);
    end
    if taxableincome>=taxbracket6 % top bracket
        tax=tax+margtaxrate6*max(taxableincome-taxbracket6,0);
    end
    aftertaxincome=taxableincome-tax;
    x=a+aftertaxincome; % cash on hand before transfers (C: coh_pre_transfer)
    % the utility floor: uMin is the utility of the floor bundle when mu=1 (C: uMin=getUtility(eMin,eMin,1)), and xMin
    % the least cash-on-hand that delivers uMin at the optimal split (C getXMin, bisection on c); fixed iteration count
    uMin=eMin^(1-nu)/(1-nu)+eMin^(1-omega)/(1-omega);
    c_lo=eMin*1e-6; c_hi=2e6;
    c_mid=eMin;
    for it=1:42
        c_mid=(c_lo+c_hi)/2;
        mtot=(mu/q)^(1/omega)*c_mid^IESfrac;
        u=c_mid^(1-nu)/(1-nu)+mu*mtot^(1-omega)/(1-omega);
        if u<uMin
            c_lo=c_mid;
        else
            c_hi=c_mid;
        end
    end
    xMin=c_mid+allocfrac*c_mid^IESfrac;
    if x<xMin
        if aprime==0
            F=uMin;
        end
    else
        c=x-aprime-q*d;
        if c>0
            F=c^(1-nu)/(1-nu)+mu*d^(1-omega)/(1-omega);
        end
    end
end

%% Death (h=2): no bequest motive; assets must be zero from then on
if h==2
    F=0;
    if aprime>0
        F=-Inf;
    end
end
%% Dead (h=3): zero
if h==3
    F=0;
end

end
