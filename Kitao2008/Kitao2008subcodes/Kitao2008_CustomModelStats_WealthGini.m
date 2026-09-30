function CustomStats=Kitao2008_CustomModelStats_WealthGini(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions)

AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Parameters,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% We just want the Gini coefficient of the asset distribution
CustomStats.GiniWealth=AllStats.A.Gini;

end