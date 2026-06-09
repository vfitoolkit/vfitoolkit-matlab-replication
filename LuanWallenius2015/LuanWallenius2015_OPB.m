function opb=LuanWallenius2015_OPB(swa,dia,agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,collar,agej_65,regime)
% Occupational pension benefit. Switches on regime:
%   regime in {0, 1} (pre-reform or regular-only reform): Swedish pre-reform
%     formula (LW2015 p.130):
%       collar = 0 blue-collar: 0.1 * (AP_55-59_best3 + 1) * BA, no claim
%                  before 65, late-claim bonus 0.7%/month up to 70.
%       collar = 1 white-collar: piecewise on wage at retirement, early-claim
%                  reduction 0.6%/month, claimable from 55 if stopped.
%       Year-factor proration applied in both.
%   regime == 2 (full reform; LW2015 p.133): DC scheme for both collars:
%       contribution per working year = 0.045 * min(earnings, 7.5) + 0.30 * max(0, earnings - 7.5)
%       benefit = capital_OPB / annuity_factor(opb_claim_age)
%       where annuity_factor is remaining model life from claim age.
%       Claim-age logic is kept the same as pre-reform; the DC annuity already
%       encodes actuarial fairness, so no extra adjustment.
% DI claimants get OPB_DI top-up instead (handled in LuanWallenius2015_OPB_DI.m),
% regardless of regime.

NEVER=0;
opb=0;

if dia>NEVER
    return;
end

% --- OPB claim age (same logic in all regimes) ---
agej_55=55-agejshifter;
if swa==NEVER
    opb_claim_agej=agej_65;
elseif collar==1
    opb_claim_agej=max(swa,agej_55);
else
    opb_claim_agej=agej_65;
end
if agej<opb_claim_agej
    return;
end

if regime<2
    % ====== Pre-reform formula (regime 0 or 1) ======

    if swa>NEVER
        yrs_worked=swa-1;
    else
        yrs_worked=agej-1;
    end
    year_factor=min(yrs_worked/30,1);

    if collar==0
        % Blue-collar: average pension points from real ages 55..59, capped at swa-1
        eff_stop_real_minus_1=(swa>NEVER)*(swa-1+agejshifter)+(swa==NEVER)*(agej-1+agejshifter);
        pp_5559_cumulator=0;
        for real_ages_5559=55:59
            hump_5559=max(0,1-((real_ages_5559-50)/30).^2);
            w_5559=w_low_i+w_peak_premium_i*hump_5559;
            pp_5559=max(0,min(w_5559*lbar,ppmax+1)-1);
            include_mask=(real_ages_5559<=eff_stop_real_minus_1);
            pp_5559_cumulator=pp_5559_cumulator+pp_5559.*include_mask;
        end
        avg_atp_5559=pp_5559_cumulator/3;
        base=0.1*(avg_atp_5559+1)*year_factor;
    else
        % White-collar: piecewise on wage at last working year
        if swa>NEVER
            j_retire=swa-1;
        else
            j_retire=agej-1;
        end
        if j_retire>=1
            hump_retire=max(0,1-((j_retire+agejshifter-50)/30)^2);
            w_retire=(w_low_i+w_peak_premium_i*hump_retire)*lbar;
        else
            w_retire=0;
        end
        bracket1=min(w_retire,7.5);
        bracket2=max(0,min(w_retire,20)-7.5);
        bracket3=max(0,min(w_retire,30)-20);
        base=(0.10*bracket1+0.65*bracket2+0.325*bracket3)*year_factor;
    end

    % Actuarial adjustment
    real_claim_age=opb_claim_agej+agejshifter;
    if collar==1
        if real_claim_age<65
            months_early=12*(65-real_claim_age);
            adj=max(1-0.006*months_early,0);
        else
            adj=1;
        end
    else
        real_capped=min(real_claim_age,70);
        if real_capped>65
            months_late=12*(real_capped-65);
            adj=1+0.007*months_late;
        else
            adj=1;
        end
    end

    opb=adj*base;
else
    % ====== Post-reform DC formula (regime 2; LW2015 p.133) ======
    % Capital accrues at 4.5% on capped earnings to 7.5BA + 30% above
    if swa>NEVER
        accrual_end=swa-1;
    else
        accrual_end=agej-1;
    end
    capital_OPB=0;
    for j=1:accrual_end
        hump_j=max(0,1-((j+agejshifter-50)/30)^2);
        w_j=w_low_i+w_peak_premium_i*hump_j;
        earnings_j=w_j*lbar;
        bracket1=min(earnings_j,7.5);
        bracket2=max(0,earnings_j-7.5);
        capital_OPB=capital_OPB+0.045*bracket1+0.30*bracket2;
    end
    annuity_factor=N_j-opb_claim_agej+1;
    opb=capital_OPB/annuity_factor;
end

end
