function CustomStats=CDGRR2003_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions)
% Calculate various calibration targets 

Params.w=(1-Params.theta)*(((Params.r+Params.delta)/(Params.theta))^(Params.theta/(Params.theta-1)));

%% First, do AllStats
% Most of the FnToEvaluate are already set up, just add a few extras
% FnsToEvaluate.K = @(l,kprime,k,s) k; %K
FnsToEvaluate.I = @(l,kprime,k,s,delta) kprime-k*(1-delta); %I
% FnsToEvaluate.L = @(l,kprime,k,s,e1,e2,e3,e4) l*(e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4)); % Efficiency hours worked: L
FnsToEvaluate.H = @(l,kprime,k,s) l; %H
% FnsToEvaluate.IncomeTaxRevenue = @(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_IncomeTaxRevenueFn(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);
% FnsToEvaluate.Pensions = @(l,kprime,k,s,J,omega) omega*(s>J); % If you are retired you earn pension omega (otherwise it is zero).
% FnsToEvaluate.EstateTaxRevenue  = @(l,kprime,k,s,J,p_gg,zlowerbar,tauE) (s>J)*(1-p_gg)*tauE*max(kprime-zlowerbar,0); % If you are retired: the probability of dying times the estate tax you would pay
FnsToEvaluate.Consumption = @(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3) CDGRR2003_ConsumptionFn(l,kprime,k,s,J,r,theta,delta,omega,e1,e2,e3,e4,a0,a1,a2,a3);

% FnsToEvaluate.Earnings = @(h,kprime,k,s,w,e1,e2,e3,e4) w*h*(e1*(s==1)+e2*(s==2)+e3*(s==3)+e4*(s==4)); 
% FnsToEvaluate.Wealth = @(h,kprime,k,s) k; % This duplicates K, but makes things easier to read as I can use 'K' in general eqm eqns and 'Wealth' in calibration targets [runtime loss is minor]

simoptions.npoints=100; % number of points for Lorenz Curve
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

%% Numerous targets are just simple ratios
Y=(AllStats.K.Mean^Params.theta)*(AllStats.L.Mean^(1-Params.theta));

CustomStats.CapitalOutputRatio=AllStats.K.Mean/Y; % K/Y
CustomStats.GovExpenditureToOutputRatio=(AllStats.IncomeTaxRevenue.Mean+AllStats.EstateTaxRevenue.Mean-AllStats.Pensions.Mean)/Y;
CustomStats.TransfersToOutputRatio=AllStats.Pensions.Mean/Y;

CustomStats.ShareOfDisposableTimeAllocatedToMarket=AllStats.H.Mean/Params.elle; % h % Working one-third of the time...???

CustomStats.zlowerbarMinus10timesAverageIncome=Params.zlowerbar-10*Y; % Not 100 percent sure this is how CDGRR2003 thought of 'average income' in this calibration target.
CustomStats.EstateTaxRevenueAsFractionOfGDP=AllStats.EstateTaxRevenue.Mean/Y;
CustomStats.EffectiveTaxRateOnAverageHHIncome=AllStats.IncomeTaxRevenue.Mean/Y; % Not clear from CDGRR2003 if they focus on income tax or all taxes. I focus on income tax.

%% Coefficient of Variance is the ratio of Standard Deviation to Mean
CustomStats.RatioOfCoeffOfVarForConsumptionToCoeffOfVarForHoursWorked=gather((AllStats.Consumption.StdDeviation/AllStats.Consumption.Mean)/(AllStats.H.StdDeviation/AllStats.H.Mean)); % Coefficient of Variation=std deviation divided by mean. 

%% Inequality statistics
% Calculate Distributions of Earnings and Wealth
CustomStats.EarningsQuintileSharesAsFraction=AllStats.Earnings.LorenzCurve([20,40,60,80,100])-AllStats.Earnings.LorenzCurve([1,21,41,61,81]);
CustomStats.EarningsTopSharesAsFraction=AllStats.Earnings.LorenzCurve([95,99,100])-AllStats.Earnings.LorenzCurve([90,95,99]);
CustomStats.WealthQuintileSharesAsFraction=AllStats.Wealth.LorenzCurve([20,40,60,80,100])-AllStats.Wealth.LorenzCurve([1,21,41,61,81]);
CustomStats.WealthTopSharesAsFraction=AllStats.Wealth.LorenzCurve([95,99,100])-AllStats.Wealth.LorenzCurve([90,95,99]);

%% Almost all the run time is to compute the two unusual model stats
% The following two model moments are somewhat unusual so required custom functions rather than using more standard VFI toolkit commands.

% Ratio of earnings of 40 year olds to 20 year olds. This is quite complicated to calculate and so required a dedicated script.
CustomStats.RatioOfEarningsOldtoYoung=CDGRR2003_RatioEarningsOldYoung(Params.NSims, StationaryDist, Policy, Params, simoptions,n_d,n_a,n_z, d_grid, a_grid,z_grid, pi_z);
% Intergenerational correlation coefficient. This is quite complicated to calculate and so required a dedicated script.
CustomStats.CrossSectionalCorrelationOfIncomeBetweenFathersAndSons=CDGRR2003_IntergenerationalEarnings(Params.NSims,StationaryDist, Policy, Params, simoptions, n_d,n_a,n_z,d_grid, a_grid,z_grid,pi_z);


end