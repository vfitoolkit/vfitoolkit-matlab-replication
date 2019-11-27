function ModelTargets=CDGRR2003_ModelTargetsFn(Params, n_d,n_a,n_z,n_p,a_grid,ReturnFn,ReturnFnParamNames, DiscountFactorParamNames,Case2_Type,PhiaprimeParamNames,FnsToEvaluateParamNames,FnsToEvaluate,GeneralEqmEqnParamNames,GEPriceParamNames,GeneralEqmEqns, vfoptions,simoptions)

Params.w=(1-Params.theta)*(((Params.r+Params.delta)/(Params.theta))^(Params.theta/(Params.theta-1)));

%% The following are just directly copied from CastanedaDiazGimenezRiosRull2003
% Since this is just a copy I have removed most of the comments here.

[e_grid,Gamma,gammastar,gammastarfull]=CastanedaDiazGimenezRiosRull2003_Create_Exog_Shock(Params);
l_grid=linspace(0,Params.elle,n_d(1))';

% Bring model into the notational conventions used by the toolkit
d_grid=[l_grid; a_grid]; % Is a 'Case 2' value function problem
z_grid=linspace(1,2*Params.J,2*Params.J)'; %(age (& determines retirement))
pi_z=Gamma;

% For 'bad' calibrations it is possible that some elements of pi_z are
% negative, if this occours then simply return that the model targets are
% all equal to zero (which will hopefully be judged as a bad calibration.
if min(min(pi_z))<0
    fprintf('Bad parameters as min(min(pi_z))<0: returning ModelTargets as all 10^5 \n')
    ModelTargets.GE_InterestRate=10^5;
    ModelTargets.GE_GovBudgetBalance=10^5;
    ModelTargets.CapitalOutputRatio=10^5;
    ModelTargets.GovExpenditureToOutputRatio=10^5;
    ModelTargets.TransfersToOutputRatio=10^5;
    ModelTargets.ShareOfDisposableTimeAllocatedToMarket=10^5;
    ModelTargets.zlowerbarMinus10timesAverageIncome=10^5;
    ModelTargets.EstateTaxRevenueAsFractionOfGDP=10^5;
    ModelTargets.EffectiveTaxRateOnAverageHHIncome=10^5;
    ModelTargets.RatioOfCoeffOfVarForConsumptionToCoeffOfVarForHoursWorked=10^5;
    ModelTargets.EarningsGini=10^5;
    ModelTargets.EarningsQuintileSharesAsFraction=[0,0,0,0,0];
    ModelTargets.EarningsTopSharesAsFraction=[0,0,0];
    ModelTargets.WealthGini=10^5;
    ModelTargets.WealthQuintileSharesAsFraction=[0,0,0,0,0];
    ModelTargets.WealthTopSharesAsFraction=[0,0,0];
    ModelTargets.RatioOfEarningsOldtoYoung=10^5;
    ModelTargets.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=10^5;
    return
end
pi_z
sum(pi_z,2)

l_p=length(n_p);
p=nan(length(GEPriceParamNames),1);
for ii=1:l_p
    p(ii)=Params.(GEPriceParamNames{ii});
end

Phi_aprimeMatrix=CastanedaDiazGimenezRiosRull2003_PhiaprimeMatrix(n_d,n_z,a_grid,Params.J,Params.zlowerbar,Params.tauE);

V0=zeros(n_a,n_z,'gpuArray');
% [p_eqm,p_eqm_index,MarketClearance]=HeteroAgentStationaryEqm_Case2(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid,Phi_aprimeMatrix, Case2_Type, ReturnFn, FnsToEvaluateFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, PriceParamNames,heteroagentoptions, simoptions, vfoptions);

fprintf('TargetsFn: Solve value fn \n')
[~, Policy]=ValueFnIter_Case2(V0, n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, Phi_aprimeMatrix, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, PhiaprimeParamNames, vfoptions);
% fprintf('TargetsFn: Some parts of the Policy function: \n')
% [Policy(2,1:10,1);...
%  Policy(2,end-9:end,1);...
%  Policy(2,1:10,4);...
%  Policy(2,end-9:end,4)]
% fprintf('TargetsFn: Value of min(min(min(Policy(2,:,:)))) =%8.2f \n', min(min(min(Policy(2,:,:)))) )
% fprintf('TargetsFn: Value of max(max(max(Policy(2,:,:)))) =%8.2f \n', max(max(max(Policy(2,:,:)))) )
% fprintf('TargetsFn: Value of min(min(min(Policy(1,:,:)))) =%8.2f \n', min(min(min(Policy(1,:,:)))) )
% fprintf('TargetsFn: Value of max(max(max(Policy(1,:,:)))) =%8.2f \n', max(max(max(Policy(1,:,:)))) )


fprintf('TargetsFn: Solve stationary dist \n')
StationaryDist=StationaryDist_Case2(Policy,Phi_aprimeMatrix,Case2_Type,n_d,n_a,n_z,pi_z,simoptions);

temp=min(min(StationaryDist));
if temp<0
    fprintf('TargetsFn: Value of min(min(StationaryDist)) =%8.2f \n', temp )
    fprintf('TargetsFn: Seems that there is some kind of rounding error that sometimes gives min(min(StationaryDist))<0, so have added line to forcibly override this \n')
    StationaryDist(StationaryDist<0)=0;
    StationaryDist=StationaryDist./sum(sum(StationaryDist));
    fprintf('TargetsFn: (Corrected) Value of min(min(StationaryDist)) =%8.2f \n', min(min(StationaryDist)) )
end
fprintf('TargetsFn: Total mass of stationary dist=%8.2f \n', sum(sum(StationaryDist)))


AggVars=EvalFnOnAgentDist_AggVars_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,2); % The 2 is for Parallel (use GPU)

% use of real() is a hack that could disguise errors, but I couldn't find why matlab was treating output as complex
GeneralEqmConditions=real(GeneralEqmConditions_Case2(AggVars,p, GeneralEqmEqns, Params,GeneralEqmEqnParamNames));

ModelTargets.GE_InterestRate=GeneralEqmConditions(1);
ModelTargets.GE_GovBudgetBalance=GeneralEqmConditions(2);

%% Calculate the Estimation Targets 

FnsToEvaluateParamNames=struct();
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_K = @(d1_val,d2_val,a_val,s_val) a_val; %K
FnsToEvaluateParamNames(2).Names={'delta'};
FnsToEvaluateFn_I = @(d1_val,d2_val,a_val,s_val,delta) d2_val-a_val*(1-delta); %I
FnsToEvaluateParamNames(3).Names={'e1','e2','e3','e4'};
FnsToEvaluateFn_L = @(d1_val,d2_val,a_val,s_val,e1,e2,e3,e4) d1_val*(e1*(s_val==1)+e2*(s_val==2)+e3*(s_val==3)+e4*(s_val==4)); % Efficiency hours worked: L
FnsToEvaluateParamNames(4).Names={};
FnsToEvaluateFn_H = @(d1_val,d2_val,a_val,s_val) d1_val; %H
FnsToEvaluateParamNames(5).Names={'J','r','theta','delta','omega','e1','e2','e3','e4','a0','a1','a2','a3'};
FnsToEvaluateFn_IncomeTaxRevenue = @(d1_val,d2_val,a_val,s_val,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_IncomeTaxRevenueFn(d1_val,d2_val,a_val,s_val,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
FnsToEvaluateParamNames(6).Names={'J','omega'};
FnsToEvaluateFn_Pensions = @(d1_val,d2_val,a_val,s_val,J,omega) omega*(s_val>J); % If you are retired you earn pension omega (otherwise it is zero).
FnsToEvaluateParamNames(7).Names={'J','p_gg','zlowerbar','tauE'};
FnsToEvaluateFn_EstateTaxRev  = @(d1_val,d2_val,a_val,s_val,J,p_gg,zlowerbar,tauE) (s_val>J)*(1-p_gg)*tauE*max(d2_val-zlowerbar,0); % If you are retired: the probability of dying times the estate tax you would pay
FnsToEvaluateParamNames(8).Names={'J','r','theta','delta','omega','e1','e2','e3','e4','a0','a1','a2','a3'};
FnsToEvaluateFn_Consumption = @(d1_val,d2_val,a_val,s_val,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_ConsumptionFn(d1_val,d2_val,a_val,s_val,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
FnsToEvaluate={FnsToEvaluateFn_K,FnsToEvaluateFn_I,FnsToEvaluateFn_L,FnsToEvaluateFn_H,FnsToEvaluateFn_IncomeTaxRevenue,FnsToEvaluateFn_Pensions,FnsToEvaluateFn_EstateTaxRev,FnsToEvaluateFn_Consumption};
AggVars=EvalFnOnAgentDist_AggVars_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,2); % The 2 is for Parallel (use GPU)

Y=(AggVars(1)^Params.theta)*(AggVars(3)^(1-Params.theta));

ModelTargets.CapitalOutputRatio=AggVars(1)/Y; % K/Y
ModelTargets.GovExpenditureToOutputRatio=(AggVars(5)+AggVars(7)-AggVars(6))/Y;
ModelTargets.TransfersToOutputRatio=AggVars(6)/Y;

ModelTargets.ShareOfDisposableTimeAllocatedToMarket=AggVars(4)/Params.elle; % h % Working one-third of the time...???

ModelTargets.zlowerbarMinus10timesAverageIncome=Params.zlowerbar-10*Y; % Not 100 percent sure this is how CDGRR2003 thought of 'average income' in this calibration target.
ModelTargets.EstateTaxRevenueAsFractionOfGDP=AggVars(7)/Y;
ModelTargets.EffectiveTaxRateOnAverageHHIncome=AggVars(5)/Y; % Not clear from CDGRR2003 if they focus on income tax or all taxes. I focus on income tax.

fprintf('Target Fn: Aggregate Variables are \n')
AggVars

FnsToEvaluateParamNames(1).Names={}; % FnsToEvaluateFn_H
FnsToEvaluateParamNames(2).Names={'J','r','theta','delta','omega','e1','e2','e3','e4','a0','a1','a2','a3'}; % FnsToEvaluateFn_Consumption
FnsToEvaluate={FnsToEvaluateFn_H, FnsToEvaluateFn_Consumption};
StationaryDist_MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,2); % The 2 is for Parallel (use GPU)

ModelTargets.RatioOfCoeffOfVarForConsumptionToCoeffOfVarForHoursWorked=(StationaryDist_MeanMedianStdDev(2,1)/StationaryDist_MeanMedianStdDev(2,3))/(StationaryDist_MeanMedianStdDev(1,1)/StationaryDist_MeanMedianStdDev(1,3)); % Coefficient of Variation=std deviation divided by mean. 

% Lorenz Curves
FnsToEvaluateParamNames(1).Names={'e1','e2','e3','e4'}; % L
FnsToEvaluateParamNames(2).Names={}; % K
FnsToEvaluateParamNames(3).Names={'J','r','theta','delta','omega','e1','e2','e3','e4','a0','a1','a2','a3'}; % FnsToEvaluateFn_Consumption
FnsToEvaluate={FnsToEvaluateFn_L,FnsToEvaluateFn_K ,FnsToEvaluateFn_Consumption}; % Note: Since we are looking at Lorenz curve of earnings we can ignore 'w' as a multiplicative scalar so will have no effect on Lorenz curve of earnings (beyond influence on d1)
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case2(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,2);
% Calculate Distributions of Earnings and Wealth
ModelTargets.EarningsGini=Gini_from_LorenzCurve(LorenzCurves(1,:));
ModelTargets.EarningsQuintileSharesAsFraction=LorenzCurves(1,[20,40,60,80,100])-LorenzCurves(1,[1,21,41,61,81]);
ModelTargets.EarningsTopSharesAsFraction=LorenzCurves(1,[95,99,100])-LorenzCurves(1,[90,95,99]);
ModelTargets.WealthGini=Gini_from_LorenzCurve(LorenzCurves(2,:));
ModelTargets.WealthQuintileSharesAsFraction=LorenzCurves(2,[20,40,60,80,100])-LorenzCurves(2,[1,21,41,61,81]);
ModelTargets.WealthTopSharesAsFraction=LorenzCurves(2,[95,99,100])-LorenzCurves(2,[90,95,99]);

% The following two model moments are somewhat unusual so required custom functions rather than using more standard VFI toolkit commands.
NSimulations=10^6;
e=[Params.e1,Params.e2,Params.e3,Params.e4,0,0,0,0];
% Ratio of earnings of 40 year olds to 20 year olds. This is quite complicated to calculate and so required a dedicated script.
ModelTargets.RatioOfEarningsOldtoYoung=CDGRR2003_RatioEarningsOldYoung(NSimulations, StationaryDist, Policy, Phi_aprimeMatrix, n_d,n_a,n_z, d_grid, Gamma, e,Params.w,Params.J);
% Intergenerational correlation coefficient. This is quite complicated to calculate and so required a dedicated script.
ModelTargets.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=CDGRR2003_IntergenerationalEarnings(NSimulations,StationaryDist, Policy, Phi_aprimeMatrix, n_d,n_a,n_z,d_grid, Gamma, e,Params.w,Params.J);

% Print out some info on what is going on:
fprintf('Targets Fn: Summary information about what is currently being evaluated \n')
Params
ModelTargets





end