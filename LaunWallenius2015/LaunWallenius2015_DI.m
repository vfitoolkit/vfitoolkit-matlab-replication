function di=LaunWallenius2015_DI(swa,dia,agej,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar,ppmax,N_j,agej_65,regime)
% Disability insurance benefit. Switches on regime:
%   regime == 0 (pre-reform; LW2015 p.131):
%       Same formula as PB but:
%         (1) no actuarial reduction for early claim
%         (2) AP computed using 3-year-pre-DI average projected forward to age 65
%         (3) benefit stays the same throughout, including after 65 auto-conversion to PB
%   regime == 1 (post-reform; LW2015 p.133):
%       DI = 0.64 * avg(w_pre_DI) * lbar
%       paid from dia onwards, including past age 65 (simplification: paper says
%       agents accrue PB capital while on DI, but we approximate that by keeping
%       the DI benefit constant -- paper p.134 notes lifetime NPV is roughly
%       unchanged across the reform for DI claimants anyway)

NEVER=0;
di=0;

if dia==NEVER || agej<dia
    return;
end

if regime==0
    % --- Pre-reform formula via AP with projection (handled inside AP when dia > 0) ---
    ap=LaunWallenius2015_AP(swa,dia,agej,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar,ppmax,N_j);

    yrs_assumed=agej_65-1;
    basic=0.96;
    earnings_supplement=0.6*ap*min(yrs_assumed/30,1);
    di=basic+earnings_supplement;
else
    % --- Post-reform: 0.64 * average wage over 3 years pre-DI ---
    proj_sum=0;
    n_proj=0;
    for j_proj=max(1,dia-3):dia-1
        real_age_pj=j_proj+agejshifter;
        w_pj=w_c2_i*real_age_pj^2+w_c1_i*real_age_pj+w_c0_i;
        proj_sum=proj_sum+w_pj;
        n_proj=n_proj+1;
    end
    if n_proj>0
        avg_pre_DI=proj_sum/n_proj;
    else
        avg_pre_DI=0;
    end
    di=0.64*avg_pre_DI*lbar;
end

end
