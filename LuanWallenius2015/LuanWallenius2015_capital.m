function capital=LuanWallenius2015_capital(accrual_end_j,agejshifter,w_low_i,w_peak_premium_i,lbar)
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
%   agejshifter, w_low_i, w_peak_premium_i, lbar: params

contrib_rate=0.185;
earnings_cap=7.5;

capital=0;
for j=1:accrual_end_j
    hump_j=max(0,1-((j+agejshifter-50)/30)^2);
    w_j=w_low_i+w_peak_premium_i*hump_j;
    capped_earnings=min(w_j*lbar,earnings_cap);
    capital=capital+contrib_rate*capped_earnings;
end

end
