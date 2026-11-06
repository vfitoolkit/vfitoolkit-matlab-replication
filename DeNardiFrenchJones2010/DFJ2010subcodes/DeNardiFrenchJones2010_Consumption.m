function c=DeNardiFrenchJones2010_Consumption(aprime,a,h,zeta,xi,r,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, normalizepsi, cfloor, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)

c=0;

if h==0 || h==1 % Alive (good or bad health)

    m_coeff=-Inf; % placeholder (just to make gpu happy)
    sigma_coeff=-Inf; % placeholder (just to make gpu happy)
    if h==0 % good health
        m_coeff=m_coeff_healthgood;
        sigma_coeff=sigma_coeff_healthgood;
    elseif h==1 % bad health
        m_coeff=m_coeff_healthbad;
        sigma_coeff=sigma_coeff_healthbad;
    end

    % Medical expenses (eqns 6,7,8 of DFJ2010)
    psi=zeta+xi; % markov + iid
    m=exp(m_coeff+sigma_coeff*(psi/normalizepsi)); % m_coeff is about 5-6, sigma_coeff is about 1-2.5
    % note: normalize psi

    % income and taxes
    taxableincome=r*a+earnings;
    % use marginal tax rates and tax brackets
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
    % after-tax income is then just
    aftertaxincome=taxableincome-tax;

    govtransfers=0;
    % Budget constraint (eqn 9 of DFJ2010)
    c=a + aftertaxincome+govtransfers-m-aprime;
    % c=a + aftertaxincome+govtransfers-aprime; % DISABLE MEDICAL EXPENSES

    % If this puts you below the consumption floor, then you can get govtransfers to reach cfloor, but only if aprime=0
    if (a+aftertaxincome-m)<cfloor && aprime==0 % [corrected 2026-09-18: was c<0, which only topped up when resources after medical expenses were NEGATIVE; DFJ2010 top up whenever they are below the floor (wealth20_Singles.c line 2257: coh=max(coh_pre_transfer,cMin))]
        % impose c=cfloor
        c=cfloor;
        % receive govtransfers to meet Consumption floor
        govtransfers=cfloor-(a+aftertaxincome-m); % note: this is just the budget constraint rearranged, but with c=cfloor and aprime=0
    end

end


end