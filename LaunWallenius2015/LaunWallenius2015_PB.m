function pb=LaunWallenius2015_PB(swa,dia,pba,agej,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar,ppmax,N_j,agej_65,regime)
% Regular pension benefit. Switches on regime:
%   regime == 0 (pre-reform Swedish PAYG-DB; LW2015 p.130):
%       PB = basic + 0.6 * AP * min(yrs/30, 1) * BA       (basic = 0.96 BA)
%       Plus actuarial adjustment for early/late claiming:
%         age 61-64 claim:  -0.5%-pts/month  (= -6%/yr)
%         age 65+ claim:    +0.7%/month     (= +8.4%/yr)
%   regime == 1 (post-reform NDC; LW2015 Section 4):
%       PB = capital_at_swa / annuity_factor(pba)
%       capital accrues at 18.5% of capped earnings each working year
%       annuity_factor(pba) = remaining model years to age 80 = N_j - pba + 1
%       Continued work after claim keeps capital growing until swa; benefit
%       is then "finalized" at capital_at_swa (LW2015 p.133).
% Paid only while pba <= agej.

NEVER=0;
pb=0;

if pba==NEVER || agej<pba
    return;
end

if regime==0
    % --- Pre-reform PAYG-DB formula ---
    ap=LaunWallenius2015_AP(swa,dia,agej,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar,ppmax,N_j);

    if swa>NEVER
        yrs_worked=swa-1;
    else
        yrs_worked=agej-1;
    end

    basic=0.96;
    earnings_supplement=0.6*ap*min(yrs_worked/30,1);
    pb_base=basic+earnings_supplement;

    real_pba=pba+agejshifter;
    if real_pba<65
        months_early=12*(65-real_pba);
        adj=1-0.005*months_early;
    elseif real_pba>65
        months_late=12*(real_pba-65);
        adj=1+0.007*months_late;
    else
        adj=1;
    end
    pb=adj*pb_base;
else
    % --- Post-reform NDC formula ---
    % accrual ends at min(swa-1, agej-1, dia-1 if on DI)
    if swa>NEVER
        accrual_end=swa-1;
    else
        accrual_end=agej-1;
    end
    if dia>NEVER
        accrual_end=min(accrual_end,dia-1);
    end
    if accrual_end<0
        accrual_end=0;
    end

    capital=LaunWallenius2015_capital(accrual_end,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar);
    annuity_factor=N_j-pba+1;
    pb=capital/annuity_factor;
end

end
