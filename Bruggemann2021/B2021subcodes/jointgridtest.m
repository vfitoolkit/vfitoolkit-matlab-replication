% Check the 

T=10;

% We need to give an initial guess for the price path
PricePath.r=Params.r*ones(T,1);
PricePath.w=Params.w*ones(T,1);
PricePath.lumpsum=Params.lumpsum*ones(T,1);

ParamPath.beta=Params.beta*ones(T,1);

n_l=10

l_grid=linspace(0,Params.maxl,n_l)'; % labor supply
% make sure lbar is a point in the grid [lbar is the fixed labor supply that entrepreneurs must provide]
[~,lbarindex]=min(abs(l_grid-Params.lbar));
l_grid(lbarindex)=Params.lbar;


%%
n_d = [n_l,n_entre];
d_grid = [l_grid; entre_grid];
tic;
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vftime_dgrid=toc
tic;
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions);
disttime_dgrid=toc
tic;
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
allstattime_dgrid=toc
tic;
PolicyValues=PolicyInd2Val_InfHorz(Policy,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
pvaltime_dgrid=toc


transpathoptions=struct();
tic;
[VPath,PolicyPath]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V, Policy, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
vpathtime=toc
tic;
AgentDistPath=AgentDistOnTransPath_InfHorz(StationaryDist, PolicyPath,n_d,n_a,n_z,pi_z,T,simoptions);
distpathtime=toc
tic;
AggVarsPath=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
aggvarpathtime=toc
tic;
PolicyValuesPath=PolicyInd2Val_InfHorz_TPath(PolicyPath,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions);
pvalpathtime=toc

%%

vfoptions.tolerance=1e-12

vfoptions.howardssparse=0
tic;
[V0,Policy0]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vftime0=toc

vfoptions.howardssparse=1
tic;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vftime1=toc

max(abs(V0(:)-V1(:)))
max(abs(Policy0(:)-Policy1(:)))
temp=Policy0(:)-Policy1(:);
sum(temp==0)
numel(temp)
sum(temp==0)/numel(temp)


%%
n_d = [n_l+1,1]; % hardcodes n_entre=2
d_grid=[[l_grid; l_grid(lbarindex)],[entre_grid(1)*ones(n_l,1); entre_grid(2)]];
tic;
[V2,Policy2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vftime=toc
tic;
StationaryDist2=StationaryDist_InfHorz(Policy2,n_d,n_a,n_z,pi_z,simoptions);
disttime=toc
tic;
AllStats2=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist2,Policy2,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
allstattime=toc
tic;
PolicyValues2=PolicyInd2Val_InfHorz(Policy2,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
pvaltime=toc

tic;
[VPath2,PolicyPath2]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V2, Policy2, Params, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionstpath);
vpathtime=toc
tic;
AgentDistPath2=AgentDistOnTransPath_InfHorz(StationaryDist2, PolicyPath2,n_d,n_a,n_z,pi_z,T,simoptions);
distpathtime=toc
tic;
AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate,AgentDistPath2,PolicyPath2,PricePath,ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);
aggvarpathtime=toc
tic;
PolicyValuesPath2=PolicyInd2Val_InfHorz_TPath(PolicyPath2,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions);
pvalpathtime=toc


%%


max(abs(V(:)-V2(:)))
max(abs(Policy(:)-Policy2(:))) % not zero
max(abs(PolicyValues(:)-PolicyValues2(:))) % but this should be zero
max(abs(StationaryDist(:)-StationaryDist2(:)))

max(abs(VPath(:)-VPath2(:)))
max(abs(PolicyPath(:)-PolicyPath2(:))) % no zero
max(abs(PolicyValuesPath(:)-PolicyValuesPath2(:))) % but this should be zero
max(abs(AgentDistPath(:)-AgentDistPath2(:)))


%%
max(abs(AllStats.Entrepreneur.Mean-AllStats2.Entrepreneur.Mean))
max(abs(AllStats.K_noncorp.Mean-AllStats2.K_noncorp.Mean))
max(abs(AllStats.A.Mean-AllStats2.A.Mean))
max(abs(AllStats.N_noncorp.Mean-AllStats2.N_noncorp.Mean))
max(abs(AllStats.N_lbar.Mean-AllStats2.N_lbar.Mean))

max(abs(AggVarsPath.Entrepreneur.Mean-AggVarsPath2.Entrepreneur.Mean))
max(abs(AggVarsPath.K_noncorp.Mean-AggVarsPath2.K_noncorp.Mean))
max(abs(AggVarsPath.A.Mean-AggVarsPath2.A.Mean))
max(abs(AggVarsPath.N_noncorp.Mean-AggVarsPath2.N_noncorp.Mean))
max(abs(AggVarsPath.N_lbar.Mean-AggVarsPath2.N_lbar.Mean))

