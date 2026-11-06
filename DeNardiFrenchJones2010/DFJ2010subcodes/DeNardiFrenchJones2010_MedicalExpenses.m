function m=DeNardiFrenchJones2010_MedicalExpenses(aprime,a,h,zeta,xi,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood,normalizepsi)

m=0;

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
    m=exp(m_coeff+sigma_coeff*(psi/normalizepsi));
    % note: normalize psi

end


end