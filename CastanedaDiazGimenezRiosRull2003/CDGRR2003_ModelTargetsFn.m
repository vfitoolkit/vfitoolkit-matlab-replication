function ModelTargets=CDGRR2003_ModelTargetsFn(Params, n_d,n_a,n_z,a_grid,ReturnFn, DiscountFactorParamNames,Case2_Type,PhiaprimeParamNames,FnsToEvaluate,GEPriceParamNames,GeneralEqmEqns, vfoptions,simoptions)

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

l_p=length(GEPriceParamNames);
p=nan(l_p,1);
for ii=1:l_p
    p(ii)=Params.(GEPriceParamNames{ii});
end

Phi_aprimeMatrix=CastanedaDiazGimenezRiosRull2003_PhiaprimeMatrix(n_d,n_z,a_grid,Params.J,Params.zlowerbar,Params.tauE);

fprintf('TargetsFn: Solve value fn \n')
[~, Policy]=ValueFnIter_Case2(n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, Phi_aprimeMatrix, Case2_Type, ReturnFn, Params, DiscountFactorParamNames, [], PhiaprimeParamNames, vfoptions);
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

AggVars=EvalFnOnAgentDist_AggVars_Case2(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);

% use of real() is a hack that could disguise errors, but I couldn't find why matlab was treating output as complex
AggVarNames=fieldnames(AggVars); % Using GeneralEqmEqns as a struct presupposes using FnsToEvaluate (and hence AggVars) as a stuct
for ii=1:length(AggVarNames)
    Params.(AggVarNames{ii})=AggVars.(AggVarNames{ii}).Mean;
end

GeneralEqmConditionsVec=real(GeneralEqmConditions_Case1_v2(GeneralEqmEqns, Params)); % Note, the GeneralEqmConditions_Case1_v2 actually has nothing to do with Case1 vs Case2

ModelTargets.GE_InterestRate=GeneralEqmConditionsVec(1);
ModelTargets.GE_GovBudgetBalance=GeneralEqmConditionsVec(2);

%% Calculate the Estimation Targets 

FnsToEvaluate.K = @(l,kprime,k,s) k; %K
FnsToEvaluate.I = @(l,kprime,k,s,delta) kprime-k*(1-delta); %I
FnsToEvaluate.L = @(l,kprime,k,s,e1,e2,e3,e4) l*(e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4)); % Efficiency hours worked: L
FnsToEvaluate.H = @(l,kprime,k,s) l; %H
FnsToEvaluate.IncomeTaxRevenue = @(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_IncomeTaxRevenueFn(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
FnsToEvaluate.Pensions = @(l,kprime,k,s,J,omega) omega*(s>J); % If you are retired you earn pension omega (otherwise it is zero).
FnsToEvaluate.EstateTaxRevenue  = @(l,kprime,k,s,J,p_gg,zlowerbar,tauE) (s>J)*(1-p_gg)*tauE*max(kprime-zlowerbar,0); % If you are retired: the probability of dying times the estate tax you would pay
FnsToEvaluate.Consumption = @(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_ConsumptionFn(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
AggVars=EvalFnOnAgentDist_AggVars_Case2(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);

Y=(AggVars.K.Mean^Params.theta)*(AggVars.L.Mean^(1-Params.theta));

ModelTargets.CapitalOutputRatio=AggVars.K.Mean/Y; % K/Y
ModelTargets.GovExpenditureToOutputRatio=(AggVars.IncomeTaxRevenue.Mean+AggVars.EstateTaxRevenue.Mean-AggVars.Pensions.Mean)/Y;
ModelTargets.TransfersToOutputRatio=AggVars.Pensions.Mean/Y;

ModelTargets.ShareOfDisposableTimeAllocatedToMarket=AggVars.H.Mean/Params.elle; % h % Working one-third of the time...???

ModelTargets.zlowerbarMinus10timesAverageIncome=Params.zlowerbar-10*Y; % Not 100 percent sure this is how CDGRR2003 thought of 'average income' in this calibration target.
ModelTargets.EstateTaxRevenueAsFractionOfGDP=AggVars.EstateTaxRevenue.Mean/Y;
ModelTargets.EffectiveTaxRateOnAverageHHIncome=AggVars.IncomeTaxRevenue.Mean/Y; % Not clear from CDGRR2003 if they focus on income tax or all taxes. I focus on income tax.

AggVarNames=fieldnames(AggVars);
fprintf('Target Fn: Aggregate Variables are \n')
for ii=1:length(AggVarNames)
    fprintf('	%s: %8.4f \n',AggVarNames{ii},AggVars.(AggVarNames{ii}).Mean)
end

FnsToEvaluate2.H=FnsToEvaluate.H;
FnsToEvaluate2.Consumption=FnsToEvaluate.Consumption;
MeanMedianStdDev=EvalFnOnAgentDist_MeanMedianStdDev_Case2(StationaryDist, Policy, FnsToEvaluate2, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);

ModelTargets.RatioOfCoeffOfVarForConsumptionToCoeffOfVarForHoursWorked=gather((MeanMedianStdDev.Consumption.StdDev/MeanMedianStdDev.Consumption.Mean)/(MeanMedianStdDev.H.StdDev/MeanMedianStdDev.H.Mean)); % Coefficient of Variation=std deviation divided by mean. 

% Lorenz Curves
FnsToEvaluate3.L=FnsToEvaluate.L;
FnsToEvaluate3.K=FnsToEvaluate.K;
LorenzCurves=EvalFnOnAgentDist_LorenzCurve_Case2(StationaryDist, Policy, FnsToEvaluate3, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid);
% Calculate Distributions of Earnings and Wealth
ModelTargets.EarningsGini=Gini_from_LorenzCurve(LorenzCurves.L);
ModelTargets.EarningsQuintileSharesAsFraction=LorenzCurves.L([20,40,60,80,100])-LorenzCurves.L([1,21,41,61,81]);
ModelTargets.EarningsTopSharesAsFraction=LorenzCurves.L([95,99,100])-LorenzCurves.L([90,95,99]);
ModelTargets.WealthGini=Gini_from_LorenzCurve(LorenzCurves.K);
ModelTargets.WealthQuintileSharesAsFraction=LorenzCurves.K([20,40,60,80,100])-LorenzCurves.K([1,21,41,61,81]);
ModelTargets.WealthTopSharesAsFraction=LorenzCurves.K([95,99,100])-LorenzCurves.K([90,95,99]);

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