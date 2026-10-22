function CustomStats=HuggettVenturaYaron2011_QuantileRatios(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,caliboptions,vfoptions,simoptions)
% CustomModelStats for the HVY2011 Section V.C sensitivity on the moments to match: the 90-10 and 75-25 quantile ratios of
% earnings by age, for the working ages 23-60 (model periods 1 to Jr-1).
% The quantiles are from the life-cycle profiles with 20 quantiles (QuantileCutoffs has the minimum, the 5,10,...,95 percent
% cutoffs, and the maximum, so the 10th, 25th, 75th and 90th percentiles are rows 3, 6, 16 and 19).

FnsToEvaluateQ.earnings=@(l,aprime,a,h,agej,Jr,g) ((1+g)^(agej-1))*h*l*(agej<Jr);
simoptionsQ=simoptions;
simoptionsQ.nquantiles=20;
AgeConditionalStatsQ=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluateQ,Parameters,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptionsQ);
QC=AgeConditionalStatsQ.earnings.QuantileCutoffs; % [21,N_j]

workages=1:1:(Parameters.Jr-1);
CustomStats.p9010=gather(QC(19,workages)./QC(3,workages));
CustomStats.p7525=gather(QC(16,workages)./QC(6,workages));

end
