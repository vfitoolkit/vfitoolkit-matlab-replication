function capital=LaunWallenius2015_capital(accrual_end_j,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar)
% Total pension capital accumulated through age accrual_end_j (inclusive),
% under the post-reform NDC contribution rule (LW2015 Section 4):
%   contribution rate = 0.185 of capped earnings each working year
%   earnings cap     = 7.5 BA per year
%
% Pure scalar implementation matching the OPB / AP scalar-loop style.
% Closed-form wage inlined.
%
% Inputs (all scalar):
%   accrual_end_j: last model period of accrual (e.g. swa-1 for stopped; agej-1 for still-working)
%   agejshifter, w_c0_i, w_c1_i, w_c2_i, lbar: params

contrib_rate=0.185;
earnings_cap=7.5;

capital=0;
for j=1:accrual_end_j
    real_age_j=j+agejshifter;
    w_j=w_c2_i*real_age_j^2+w_c1_i*real_age_j+w_c0_i;
    capped_earnings=min(w_j*lbar,earnings_cap);
    capital=capital+contrib_rate*capped_earnings;
end

end
