function net_surplus=LuanWallenius2015_BudgetFn(d,kprime,k,hist_idx,h,agej,r,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_65,regime)
% Net per-agent contribution to the government budget at current state:
%   net_surplus = (income_tax + ss_tax + consumption_tax) - (PB + DI + OPB + OPB_DI)
%
% At general-equilibrium budget balance, mean(net_surplus) = T (the lump-sum
% transfer per agent), because gov receipts = gov spending including T*mass.
% So the T-clearing condition is mean(net_surplus) = T, which is what the GE
% loop iterates on.
%
% Mirrors the consumption / tax / benefit logic from LuanWallenius2015_ReturnFn,
% with the wage formula inlined the same way.

NEVER=0;
n_pba_options=21;
n_no_di=57*n_pba_options;       % 1197

% --- Decode hist_idx -> (swa, dia, pba) ---
if hist_idx<=n_no_di
    swa=floor((hist_idx-1)/n_pba_options);
    pba_idx=mod(hist_idx-1,n_pba_options)+1;
    if pba_idx==1
        pba=NEVER;
    else
        pba=pba_idx+35;
    end
    dia=NEVER;
else
    dia=hist_idx-n_no_di;
    swa=dia;
    pba=NEVER;
end

% --- Earnings (only if working this period) ---
if d==1 || d==5
    hump=max(0,1-((agej+agejshifter-50)/30)^2);
    w_age=w_low_i+w_peak_premium_i*hump;
    earnings=w_age*lbar;
else
    earnings=0;
end

% --- Benefits ---
PB    =LuanWallenius2015_PB    (swa,dia,pba,agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,agej_65,regime);
DI    =LuanWallenius2015_DI    (swa,dia,    agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,agej_65,regime);
OPB   =LuanWallenius2015_OPB   (swa,dia,    agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,collar_i,agej_65,regime);
OPB_DI=LuanWallenius2015_OPB_DI(    dia,        agejshifter,w_low_i,w_peak_premium_i,lbar,                collar_i);
total_benefits=PB+DI+OPB+OPB_DI;

% --- Taxes ---
Y=earnings+PB+DI+OPB+OPB_DI;
if Y<=tax_thresh
    income_tax=tau_l1*Y;
else
    income_tax=tau_l1*tax_thresh+tau_l2*(Y-tax_thresh);
end
ss_tax=tau_ss*earnings;

% --- Consumption (from budget; depends on T) ---
after_tax=Y-income_tax-ss_tax;
gross_expend=(1+r)*k+after_tax+T-kprime;
if gross_expend<=0
    c=0;
else
    c=gross_expend/(1+tau_c);
end
consumption_tax=tau_c*c;

% --- Net surplus per agent ---
total_tax_revenue=income_tax+ss_tax+consumption_tax;
net_surplus=total_tax_revenue-total_benefits;

end
