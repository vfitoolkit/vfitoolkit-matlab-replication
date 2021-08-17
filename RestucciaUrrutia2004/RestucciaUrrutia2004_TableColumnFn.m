function TableColumn=RestucciaUrrutia2004_TableColumnFn(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions)
% Creates all the output necessary for one of the columns of Tables 5, 6, 7, or 8

% Calculate the General Eqm for the current parameter values
[p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case2_FHorz(jequaloneDist,AgeWeightParamNames,n_d, n_a, n_z, N_j, n_p, AgeDependentGridParamNames, d_gridfn, a_gridfn, z_gridfn,Phi_aprime, Case2_Type, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
Params.g=p_eqm.g;

[V, Policy]=ValueFnIter_Case2_FHorz(n_d,n_a,n_z,N_j, d_gridfn, a_gridfn, z_gridfn, AgeDependentGridParamNames, Phi_aprime, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, vfoptions);
StationaryDist=StationaryDist_FHorz_Case2(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,d_gridfn, a_gridfn, z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,Params,PhiaprimeParamNames,simoptions);

%% First six rows are just those that were needed for Table 2
% TableColumn.IntergenerationalCorrelationOfInnateAbility
% TableColumn.IntergenerationalCorrelationOfAcquiredAbility
% TableColumn.IntergenerationalCorrelationOfEarnings
% TableColumn.StdDevLogInnateAbility
% TableColumn.StdDevLogAcquiredAbility
% TableColumn.StdDevLogEarnings

% For Table 8 we also need intergenerational correlations of 'educational
% attainment (college completion)' and consumption
% TableColumn.IntergenerationalCorrelationOfEducationalAttainment
% TableColumn.IntergenerationalCorrelationOfConsumption

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'agej','w'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w) log(w*a1_val); % Earnings
FnsToEvaluateParamNames(2).Names={'agej'};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) log(a2_val); % bhat, (log) aquired ability
FnsToEvaluateParamNames(3).Names={'agej'};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) log(z1_val); % b, (log) innate ability
FnsToEvaluateParamNames(4).Names={'agej'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) agej; % age
FnsToEvaluateParamNames(5).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_5 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % indicator for college completion
FnsToEvaluateParamNames(6).Names={'agej','w','g','gamma','tau','kappa0','kappa1','plowerbar', 'pupperbar','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_6 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,g,gamma,tau,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_Consumption(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,g,gamma,tau,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1); % consumption
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6};

% First, the standard deviations in the cross section. Since the number for
% earnings is same as that in Table 1 this has presumably also been
% calcuated as conditional on agej==2. I assume this is true of all the
% numbers here.
simoptions.parallel=3;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
simoptions.parallel=2;
% The following are elements of Table 2 that are really age-contional moments
TableColumn.StdDevLogEarnings=gather(AgeConditionalStats(1).StdDeviation(2));
TableColumn.StdDevLogAcquiredAbility=gather(AgeConditionalStats(2).StdDeviation(2));
TableColumn.StdDevLogInnateAbility=gather(AgeConditionalStats(3).StdDeviation(2));

% Second, the intergenerational correlations.
PanelForCalculatingIntergenCorr=SimPanelValues_FHorz_Case2(StationaryDist,Policy, FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,PhiaprimeParamNames, simoptions);
% Earnings
X1=PanelForCalculatingIntergenCorr(1,1:end-1,:).*(PanelForCalculatingIntergenCorr(4,1:end-1,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(1,2:end,:).*(PanelForCalculatingIntergenCorr(4,2:end,:)==1);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
TableColumn.IntergenerationalCorrelationOfEarnings=regcoeff(2);
% Acquired Ability (only measured when j==2, so have to do slightly differently)
X1=PanelForCalculatingIntergenCorr(2,1:end-2,:).*(PanelForCalculatingIntergenCorr(4,1:end-2,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(2,3:end,:).*(PanelForCalculatingIntergenCorr(4,3:end,:)==2);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
TableColumn.IntergenerationalCorrelationOfAcquiredAbility=regcoeff(2); % corr(X1,Y,'rows','pairwise')
% Innate Ability
X1=PanelForCalculatingIntergenCorr(3,1:end-1,:).*(PanelForCalculatingIntergenCorr(4,1:end-1,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(3,2:end,:).*(PanelForCalculatingIntergenCorr(4,2:end,:)==1);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
TableColumn.IntergenerationalCorrelationOfInnateAbility=regcoeff(2); % corr(X1,Y,'rows','pairwise')

% Educational Attainment (college completion) (only measured when j==2, so have to do slightly differently)
X1=PanelForCalculatingIntergenCorr(5,1:end-2,:).*(PanelForCalculatingIntergenCorr(4,1:end-2,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(5,3:end,:).*(PanelForCalculatingIntergenCorr(4,3:end,:)==2);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
TableColumn.IntergenerationalCorrelationOfEducationalAttainment=regcoeff(2); % corr(X1,Y,'rows','pairwise')
% Consumption(only measured when j==2??? (presumably?), so have to do slightly differently)
X1=PanelForCalculatingIntergenCorr(6,1:end-2,:).*(PanelForCalculatingIntergenCorr(4,1:end-2,:)==2);
X1=reshape(X1,[numel(X1),1]); X1(X1==0)=nan;
Y=PanelForCalculatingIntergenCorr(6,3:end,:).*(PanelForCalculatingIntergenCorr(4,3:end,:)==2);
Y=reshape(Y,[numel(Y),1]); Y(Y==0)=nan;
regcoeff=regress(Y,[ones(numel(X1),1),X1]);
TableColumn.IntergenerationalCorrelationOfConsumption=regcoeff(2); % corr(X1,Y,'rows','pairwise')



%% 
% TableColumn.PrivateEarlyEducDivGDP: e/Y
% TableColumn.AverageCollegePremium

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'agej','plowerbar', 'pupperbar','nlowerbar','nupperbar','psi0','psi1'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar, pupperbar,nlowerbar,nupperbar,psi0,psi1) RestucciaUrrutia2004_HFn(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,plowerbar, pupperbar,nlowerbar,nupperbar,psi0,psi1); % H, Human capital
FnsToEvaluateParamNames(2).Names={'agej','w','kappa0','kappa1','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_FnkappaF(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1); % kappaF (a3_val==1, is that the state s=1 for old)
FnsToEvaluateParamNames(3).Names={'agej'};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (agej==2)*(1-a3_val); % fraction of j==2 population with s==0: Table 1 (i) fraction of non-college
FnsToEvaluateParamNames(4).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % fraction college dropout
FnsToEvaluateParamNames(5).Names={'agej','gamma','g'};
FnsToEvaluateFn_5 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,gamma,g) (agej==1)*((d1_val/z1_val)^(1/gamma)-g); % e=(bhat/b)^(1/gamma)-g for age j==1, private expenditure in early education
FnsToEvaluateParamNames(6).Names={'agej','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_6 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_FnF(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,nlowerbar,nupperbar,f,psi0,psi1); % F (a3_val==1, is that the state s=1 for old)
% The next six are about getting the wages for non-college, college dropout, and college so as to calculate the (completed and dropout) college premium
FnsToEvaluateParamNames(7).Names={'agej','w'};
FnsToEvaluateFn_7 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w) w*a1_val*(a3_val==0)*(agej==2); % Wage for non-college
FnsToEvaluateParamNames(8).Names={'agej'};
FnsToEvaluateFn_8 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (a3_val==0)*(agej==2); % Fraction non-college
FnsToEvaluateParamNames(9).Names={'agej','w','psi0','psi1','plowerbar'};
FnsToEvaluateFn_9 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,psi0,psi1,plowerbar) w*plowerbar*a1_val*(a3_val==1)*(agej==2)*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % Wage for college dropout
FnsToEvaluateParamNames(10).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_10 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (a3_val==1)*(agej==2)*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % Fraction college dropout
FnsToEvaluateParamNames(11).Names={'agej','w','psi0','psi1','pupperbar'};
FnsToEvaluateFn_11 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,psi0,psi1,pupperbar) w*pupperbar*a1_val*(a3_val==1)*(agej==2)*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % Wage for college completed
FnsToEvaluateParamNames(12).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_12 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (a3_val==1)*(agej==2)*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % Fraction college completed
FnsToEvaluateParamNames(13).Names={};
FnsToEvaluateFn_13 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a1_val; % Human capital
FnsToEvaluateParamNames(14).Names={'agej','w','g','gamma','tau','kappa0','kappa1','plowerbar', 'pupperbar','nlowerbar','nupperbar','f','psi0','psi1'};
FnsToEvaluateFn_14 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,g,gamma,tau,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1) RestucciaUrrutia2004_Consumption(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,w,g,gamma,tau,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1); % consumption
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8,FnsToEvaluateFn_9,FnsToEvaluateFn_10,FnsToEvaluateFn_11,FnsToEvaluateFn_12,FnsToEvaluateFn_13,FnsToEvaluateFn_14};

MomentsForTables=gather(EvalFnOnAgentDist_AggVars_FHorz_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_gridfn, a_gridfn, z_gridfn, simoptions, AgeDependentGridParamNames));

GDP=Params.A*MomentsForTables(1); % Y=A*H from firms production function.

TableColumn.PrivateEarlyEducDivGDP=MomentsForTables(5)/GDP; % e/Y

% Average wage of non-college
AvgWage_NoCollege=MomentsForTables(7)/MomentsForTables(8);
% Average wage of college dropout
AvgWage_CollegeDropout=MomentsForTables(9)/MomentsForTables(10);
% Average wage of college
AvgWage_CollegeCompleted=MomentsForTables(11)/MomentsForTables(12);
TableColumn.AverageDropoutPremium=AvgWage_CollegeDropout/AvgWage_NoCollege;
TableColumn.AverageCollegePremium=AvgWage_CollegeCompleted/AvgWage_NoCollege;

% Couple of things for Table 8
TableColumn.ExpendAsPercentOfGDP_PrivateEarlyEducation=100*TableColumn.PrivateEarlyEducDivGDP; % e/Y % note: just renaming/duplicating
TableColumn.ExpendAsPercentOfGDP_PublicEarlyEducation=100*Params.g/GDP; % g/Y
TableColumn.ExpendAsPercentOfGDP_PrivateCollegeEducation=100*(MomentsForTables(6)-MomentsForTables(2))/GDP; % (F-kappaF)/Y (note: kappaF is not 'kappa times F')
TableColumn.ExpendAsPercentOfGDP_PublicCollegeEducation=100*MomentsForTables(2)/GDP;  % kappaF/Y
TableColumn.CollegeEnrollmentRate=1-MomentsForTables(8);
TableColumn.CollegeDropoutRate=MomentsForTables(10)/(1-MomentsForTables(8));
TableColumn.AggregateHumanCapital=MomentsForTables(13);
TableColumn.AggregateConsumption=MomentsForTables(14);


%% 
% TableColumn.CollegeEnrollmentRate

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'w'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,w) w*a1_val; % Earnings
FnsToEvaluateParamNames(2).Names={};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a2_val; % bhat, aquired ability.
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) z1_val; % b, innate ability
FnsToEvaluateParamNames(4).Names={'agej'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) agej; % age
FnsToEvaluateParamNames(5).Names={'agej','gamma','g'};
FnsToEvaluateFn_5 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,gamma,g) (agej==1)*((d1_val/z1_val)^(1/gamma)-g); % e=(bhat/b)^(1/gamma)-g for age j==1, private expenditure in early education
FnsToEvaluateParamNames(6).Names={'agej'};
FnsToEvaluateFn_6 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (a3_val==1)*(agej==2); % college enrollment decision
FnsToEvaluateParamNames(7).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_7 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % indicator for college completion
FnsToEvaluateParamNames(8).Names={};
FnsToEvaluateFn_8 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) d1_val; % bhat, aquired ability (again, but in a way that is easily paired with the young parent HH)
FnsToEvaluateParamNames(9).Names={};
FnsToEvaluateFn_9 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a1_val; % Human capital
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8,FnsToEvaluateFn_9};

simoptions.parallel=3;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
simoptions.parallel=2;

TableColumn.CollegeEnrollmentRate=AgeConditionalStats(6).Mean(2);

%% Correlation with log acquired ability (for Table 7)

simoptions.agegroupings=[1,N_j]; % Do each age seperately (note, this just evaluates to [1,2])
simoptions.crosssectioncorrelation=1;

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'w'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,w) log(1+w*a1_val); %  (log) Earnings % The "1+" is because otherwise there are some zeros which when log() is taken give -Inf and the correlations are therefore NaN
FnsToEvaluateParamNames(2).Names={};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) log(a2_val); % bhat,  (log) aquired ability.
FnsToEvaluateParamNames(3).Names={'agej'};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (a3_val==1)*(agej==2); % college enrollment decision
FnsToEvaluateParamNames(4).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % indicator for college completion
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4};
simoptions.parallel=3;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case2(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn, simoptions, AgeDependentGridParamNames);
simoptions.parallel=2;

TableColumn.Correlationwithlogacqabilt_collegeenroll=AgeConditionalStats(2,3).CrossSectionalCorrelation(2);
TableColumn.Correlationwithlogacqabilt_educattain=AgeConditionalStats(2,4).CrossSectionalCorrelation(2);
TableColumn.Correlationwithlogacqabilt_logearnings=AgeConditionalStats(2,1).CrossSectionalCorrelation(1);


%% Figures 3 and 4

% Many of these are same as for Table 2, but now without log (actually
% given we are just interested in terciles I could still use log if I
% wanted)
FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={'w'};
FnsToEvaluateFn_1 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,w) w*a1_val; % Earnings
FnsToEvaluateParamNames(2).Names={};
FnsToEvaluateFn_2 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a2_val; % bhat, aquired ability.
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluateFn_3 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) z1_val; % b, innate ability
FnsToEvaluateParamNames(4).Names={'agej'};
FnsToEvaluateFn_4 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) agej; % age
FnsToEvaluateParamNames(5).Names={};
FnsToEvaluateFn_5 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val) a1_val; % Human capital
FnsToEvaluateParamNames(6).Names={'agej'};
FnsToEvaluateFn_6 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej) (a3_val==1)*(agej==2); % college enrollment decision
FnsToEvaluateParamNames(7).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_7 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (agej==2)*a3_val*(z2_val>=min(psi0*(1+a2_val)^psi1,1)); % indicator for college completion
FnsToEvaluateParamNames(8).Names={'agej','psi0','psi1'};
FnsToEvaluateFn_8 = @(d1_val,d2_val,a1_val,a2_val,a3_val,z1_val,z2_val,z3_val,agej,psi0,psi1) (a3_val==1)*(agej==2)*(z2_val<min(psi0*(1+a2_val)^psi1,1)); % Fraction college dropout
FnsToEvaluate={FnsToEvaluateFn_1,FnsToEvaluateFn_2,FnsToEvaluateFn_3,FnsToEvaluateFn_4,FnsToEvaluateFn_5,FnsToEvaluateFn_6,FnsToEvaluateFn_7,FnsToEvaluateFn_8};

PanelForDrawingFigure=SimPanelValues_FHorz_Case2(StationaryDist,Policy, FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,N_j,d_gridfn,a_gridfn,z_gridfn,AgeDependentGridParamNames, Phi_aprime,Case2_Type,PhiaprimeParamNames, simoptions);

% For Figure
Indicator_OldHH=shiftdim((PanelForDrawingFigure(4,:,:)==2),1);
% Calculate the quintile cutoffs for old parent earnings
Earnings=shiftdim(PanelForDrawingFigure(1,:,:),1); % redundant as already calculated above
Earnings=Earnings(Indicator_OldHH);
OldEarningsTercileCutoffs=prctile(Earnings,[20,40,60,80]);
% The three variables that we wish to graph
ChildInnateAbility=shiftdim(PanelForDrawingFigure(3,:,:),1);
ChildInnateAbility=ChildInnateAbility(Indicator_OldHH);
ChildAcquiredAbility=shiftdim(PanelForDrawingFigure(2,:,:),1);
ChildAcquiredAbility=ChildAcquiredAbility(Indicator_OldHH);
HumanCapital=shiftdim(PanelForDrawingFigure(5,:,:),1);
HumanCapital=HumanCapital(Indicator_OldHH);
% Set up indicators based on the quintiles
Indicator_EarningsQuintile1=(Earnings<OldEarningsTercileCutoffs(1));
Indicator_EarningsQuintile2=(Earnings>=OldEarningsTercileCutoffs(1)).*(Earnings<OldEarningsTercileCutoffs(2));
Indicator_EarningsQuintile3=(Earnings>=OldEarningsTercileCutoffs(2)).*(Earnings<OldEarningsTercileCutoffs(3));
Indicator_EarningsQuintile4=(Earnings>=OldEarningsTercileCutoffs(3)).*(Earnings<OldEarningsTercileCutoffs(4));
Indicator_EarningsQuintile5=(Earnings>=OldEarningsTercileCutoffs(4));

% Is also off of old HH earnings quintiles. So just need to get the
% variables themselves and can reuse all the indicators from Figure 1
CollegeEnrollment=shiftdim(PanelForDrawingFigure(6,:,:),1);
CollegeEnrollment=CollegeEnrollment(Indicator_OldHH);
CollegeCompletion=shiftdim(PanelForDrawingFigure(7,:,:),1);
CollegeCompletion=CollegeCompletion(Indicator_OldHH);
CollegeDropout=shiftdim(PanelForDrawingFigure(8,:,:),1);
CollegeDropout=CollegeDropout(Indicator_OldHH);

% Note CollegeCompletion and CollegeDropout are fractions of entire
% population; the 'CollegeDropoutRate' below is dropout as fraction of
% enrolled.

TableColumn.FigureData.CollegeEnrollment=100*[mean(CollegeEnrollment(logical(Indicator_EarningsQuintile1)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile2)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile3)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile4)));...
    mean(CollegeEnrollment(logical(Indicator_EarningsQuintile5)))];

% College completion not actually used in figures
TableColumn.FigureData.CollegeCompletionRate=100*[mean(CollegeCompletion(logical(Indicator_EarningsQuintile1)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile1)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile2)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile2)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile3)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile3)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile4)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile4)));...
    mean(CollegeCompletion(logical(Indicator_EarningsQuintile5)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile5)))];

TableColumn.FigureData.CollegeDropoutRate=100*[mean(CollegeDropout(logical(Indicator_EarningsQuintile1)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile1)));...
    mean(CollegeDropout(logical(Indicator_EarningsQuintile2)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile2)));...
    mean(CollegeDropout(logical(Indicator_EarningsQuintile3)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile3)));...
    mean(CollegeDropout(logical(Indicator_EarningsQuintile4)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile4)));...
    mean(CollegeDropout(logical(Indicator_EarningsQuintile5)))./mean(CollegeEnrollment(logical(Indicator_EarningsQuintile5)))];




