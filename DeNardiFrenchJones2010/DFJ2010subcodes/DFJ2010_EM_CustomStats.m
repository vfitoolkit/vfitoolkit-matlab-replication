function CustomStats=DFJ2010_EM_CustomStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,estimoptions,vfoptions,simoptions) %#ok<INUSD>
% DFJ2010_EM_CustomStats: estimoptions.CustomModelStats for the endogenous medical-expense model, the autocorrelations of log
% out-of-pocket medical expenses in DFJ2010's OWN definition (gmmcrit51.g simcrrl; online appendix eqns A4-A5): over the
% people observed at both dates, the mean of the product of the two dates' log expenses each centred on the FULL cell's mean
% and divided by the full cell's standard deviation at that date (the toolkit's AutoCorrelation centres and scales by the
% pair population's own means and sds; the data rows ac1dfj/ac2dfj of DFJ2010_MedexTargets are this definition).
% Called by the cohort estimator once per cohort on that cohort's distribution (simoptions.jequaloneDistAge is its entry age),
% the same signature as without cohorts. The expenses are 2-year averages (simoptions.ntimeperiods=2, DFJ2010 average adjacent
% simulated years), logged (simoptions.logwhichcombos=1), by PI half (restrictions aliveh1/aliveh2, which must hold at every
% age of both windows), at horizon 2 (consecutive waves, DFJ's first autocorrelation) and horizon 4 (DFJ's second).
% From the toolkit outputs: E_pair[(x_j-m_j)(x_{j+K}-m_{j+K})] = AutoCovariance_horizonK + (PairMean_j-m_j)(PairMean_jplusk-m_{j+K}),
% with m and s the restricted full-cell Mean and StdDeviation of the log window at ages j and j+K (computed at exactly those ages).
% Returns, index = window start age (NaN where not computed):
%   AC1_h1, AC1_h2   1 x N_j-2   horizon 2, PI half 1 and 2
%   AC2_h1, AC2_h2   1 x N_j-4   horizon 4
% Only the window start ages the targets use are computed (the cohort's age in 1997, 1999, 2001 and 2003; see DFJ2010_EMTargets).
j0=simoptions.jequaloneDistAge;
Fns.oopmedex=FnsToEvaluate.oopmedex;
simoptions_a=simoptions;
simoptions_a.conditionalrestrictions=struct(); % the PI halves only: pages unrestricted, aliveh1, aliveh2
simoptions_a.conditionalrestrictions.aliveh1=simoptions.conditionalrestrictions.aliveh1;
simoptions_a.conditionalrestrictions.aliveh2=simoptions.conditionalrestrictions.aliveh2;
simoptions_a.ntimeperiods=2;
simoptions_a.logwhichcombos=1;
simoptions_a.timehorizons=[2,4];
simoptions_a.whichcombos=zeros(1,N_j,3);
origins=j0+[1,3,5,7]; origins=origins(origins<=N_j-3);
simoptions_a.whichcombos(1,origins,2:3)=1;
ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,Fns,Parameters,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_a);
CustomStats=struct();
for hh=1:2
    R=ACP.(['aliveh',num2str(hh)]).oopmedex_2periodAverage_log;
    m=gather(R.Mean); s=gather(R.StdDeviation);
    for K=[2,4]
        num=gather(R.(['AutoCovariance_horizon',num2str(K)]))+(gather(R.(['PairMean_j_horizon',num2str(K)]))-m(1:N_j-K)).*(gather(R.(['PairMean_jplusk_horizon',num2str(K)]))-m(1+K:N_j));
        CustomStats.(['AC',num2str(K/2),'_h',num2str(hh)])=num./(s(1:N_j-K).*s(1+K:N_j));
    end
end
end
