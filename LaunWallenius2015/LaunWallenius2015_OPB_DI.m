function opb_di=LaunWallenius2015_OPB_DI(dia,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar,collar)
% Occupational pension DI top-up, Swedish pre-reform (LW2015 p.131; LW2016 Tdioc.m).
%   Blue-collar (collar=0):  0.15 * min(w_preDI, 7.5)
%   White-collar (collar=1): 0.15 * min(w_preDI, 7.5)
%                          + 0.65 * max(0, min(w_preDI, 20) - 7.5)
%                          + 0.325 * max(0, min(w_preDI, 30) - 20)
% Pre-DI wage = wage at last working year (model age dia-1) times lbar.
% Paid alongside DI benefit while on DI; LW2016 code continues paying it
% after the 65-DI->PB conversion as well.
%
% Wage formula inlined (closed form; regression coefficients sent by Tobias Laun via email).

NEVER=0;
opb_di=0;
if dia==NEVER
    return;
end

if dia>=2
    j_pre_DI=dia-1;
    real_age_pre_DI=j_pre_DI+agejshifter;
    w_pre_DI=(w_c2_i*real_age_pre_DI^2+w_c1_i*real_age_pre_DI+w_c0_i)*lbar;
else
    w_pre_DI=0;
end

bracket1=min(w_pre_DI,7.5);
bracket2=max(0,min(w_pre_DI,20)-7.5);
bracket3=max(0,min(w_pre_DI,30)-20);

if collar==0
    opb_di=0.15*bracket1;
else
    opb_di=0.15*bracket1+0.65*bracket2+0.325*bracket3;
end

end
