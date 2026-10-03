function ssb=HV2000_SSBenefitFn_M4(d,aprime,a,ebar,z,e,agej,Jr,g,b_common,ybar_mean,w,haveSS)
% Model 4 wrapper: SS benefit doesn't depend on e but FnsToEvaluate
% signature must include e when n_e>0. Delegates to shared HV2000_SSBenefitFn.

ssb=HV2000_SSBenefitFn(d,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean,w,haveSS);

end
