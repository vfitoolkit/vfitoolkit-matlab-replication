function F=LuanWallenius2015_ReturnFn(d,kprime,k,hist_idx,h,agej,r,b1,b2,b3,b4,b5,tax_thresh,tau_l1,tau_l2,tau_ss,tau_c,T,w_low_i,w_peak_premium_i,collar_i,lbar,ppmax,N_j,agejshifter,agej_pb_first,agej_di_max,agej_65,regime)
% Action space: (d, kprime, k, hist_idx, h, ...)
%   d:        action code 1..5 (labor + claim decisions; see aprimeFn header)
%   kprime:   next-period assets (standard endo aprime)
%   k:        current assets
%   hist_idx: experienceasset, encodes (swa, dia, pba) -- see driver header
%   h:        exogenous health, integer 1..5 (1=very bad, 5=very good)

F=-Inf;

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

% --- Action feasibility ---
on_DI =(dia>NEVER && agej>=dia && agej<agej_65);
on_PB =(pba>NEVER && agej>=pba);

% d=1: keep working. Requires swa==NEVER, not on DI.
% d=2: stop (or stay stopped), no claim change. Always feasible if not on DI in a forced-work way (DI agents stay forced to d=2).
% d=3: apply for DI. Requires h<=2, dia==NEVER, agej<=agej_di_max, pba==NEVER, not on PB.
% d=4: claim PB and stop. Requires pba==NEVER, agej>=agej_pb_first, dia==NEVER.
% d=5: claim PB and keep working. Requires pba==NEVER, agej>=agej_pb_first, swa==NEVER, dia==NEVER.

if on_DI && d~=2
    return; % DI claimants must take d=2 (no work, no new claim)
end
if d==1
    if swa>NEVER
        return; % work absorbing
    end
elseif d==3
    if h>2 || dia>NEVER || agej>agej_di_max || pba>NEVER
        return;
    end
elseif d==4
    if pba>NEVER || agej<agej_pb_first || dia>NEVER
        return;
    end
elseif d==5
    if pba>NEVER || agej<agej_pb_first || swa>NEVER || dia>NEVER
        return;
    end
end

% --- Labor & earnings this period ---
if d==1 || d==5
    labor=lbar;
    hump=max(0,1-((agej+agejshifter-50)/30)^2);
    w_age=w_low_i+w_peak_premium_i*hump;
    earnings=w_age*labor;
else
    labor=0;
    earnings=0;
end

% --- Benefits (functions of current state, not current d) ---
PB    =LuanWallenius2015_PB    (swa,dia,pba,agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,agej_65,regime);
DI    =LuanWallenius2015_DI    (swa,dia,    agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,agej_65,regime);
OPB   =LuanWallenius2015_OPB   (swa,dia,    agej,agejshifter,w_low_i,w_peak_premium_i,lbar,ppmax,N_j,collar_i,agej_65,regime);
OPB_DI=LuanWallenius2015_OPB_DI(    dia,        agejshifter,w_low_i,w_peak_premium_i,lbar,                collar_i);

% --- Taxes (progressive labor income tax on Y; payroll tax on earnings only) ---
Y=earnings+DI+PB+OPB+OPB_DI;
if Y<=tax_thresh
    tax_y=tau_l1*Y;
else
    tax_y=tau_l1*tax_thresh+tau_l2*(Y-tax_thresh);
end
tax_ss=tau_ss*earnings;
after_tax=Y-tax_y-tax_ss;

% --- Budget: (1+tau_c)*c + kprime = (1+r)*k + after_tax + T ---
gross_expend=(1+r)*k+after_tax+T-kprime;
if gross_expend<=0
    return;
end
c=gross_expend/(1+tau_c);

% --- Utility ---
if h==1
    b=b1;
elseif h==2
    b=b2;
elseif h==3
    b=b3;
elseif h==4
    b=b4;
else
    b=b5;
end

F=log(c)-b*labor+h;

end
