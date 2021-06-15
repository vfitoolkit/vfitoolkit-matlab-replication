function output=HopenhaynRogerson1993_Fn(ImposeFootnote5, Params,n_a, n_z,a_grid,z_grid,pi_z, ReturnFn, ReturnFnParamNames, DiscountFactorParamNames, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions)
% Largely just a copy-paste of main HopenhaynRogerson1993 code, except
% towards the end when it contains all the calculations of model output
% required for Tables 3, 4 and 5.
% Only substantial difference is that here we treat 'ce' (fixed cost of
% entry) as a parameter, and 'p' as a general equilibrium price to be
% calculated.

n_d=0; % None.
d_grid=[];

FnsToEvaluateParamNames(1).Names={'alpha'};
% Note: With entry-exit the mass of the distribution of agents often
% matters. So it becomes an extra input arguement in all functions to be evaluated.
FnsToEvaluateFn_1 = @(aprime_val,a_val,z_val,agentmass,alpha) z_val*(aprime_val^alpha); % Real output
FnsToEvaluate={FnsToEvaluateFn_1};

GEPriceParamNames={'p','Ne'};
% heteroagentoptions.specialgeneqmcondn={0,'entry'};
GeneralEqmEqnParamNames(1).Names={'A'};
GeneralEqmEqn_1 = @(AggVars,GEprices,A) A/AggVars-GEprices(1); %The requirement that the price is determined by the demand eqn (or equivalently, can think of this as goods market clearance)
GeneralEqmEqnParamNames(2).Names={'ce','beta'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,ce,beta) beta*EValueFn-GEprices(1)*ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
if ImposeFootnote5==1
    GeneralEqmEqnParamNames(2).Names={'ce'};
    GeneralEqmEqn_Entry = @(EValueFn,GEprices,ce) EValueFn-GEprices(1)*ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
end
% The general eqm conditions look slightly different to more standard @(EValueFn,p,params)
% This is because 'p' is the name of a parameter, and so have used 'GEprices'
% instead of my usual 'p' to refer to the general equilibrium prices (parameter 'p' is actually GEprices(1))
GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_Entry};

n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);

Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;

[V,Policy,ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.zeta=1-ExitPolicy;
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);

% The following four are not required for anything, but will output them anyway
output.V=V;
output.Policy=Policy;
output.ExitPolicy=ExitPolicy;
output.StationaryDist=StationaryDist;

%% Now that the stationary equilibrium has been found, generate the statistics required for Tables 3, 4 and 5.

%% Table 3 stats
output.price=Params.p;

FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_Emp = @(aprime_val,a_val,z_val,AgentDistMass) aprime_val; % Employment
FnsToEvaluateParamNames(2).Names={};
FnsToEvaluateFn_Hiring = @(aprime_val,a_val,z_val,AgentDistMass) (aprime_val-a_val)*(aprime_val>a_val); % Hiring (need to add the 'firm entry' which involves hiring a single worker)
FnsToEvaluateParamNames(3).Names={};
FnsToEvaluateFn_Firing = @(aprime_val,a_val,z_val,AgentDistMass) -(aprime_val-a_val)*(aprime_val<a_val); % Firing (need to add the 'firm exits' which involve firing all remaing workers)
FnsToEvaluateParamNames(4).Names={'tau','w'};
FnsToEvaluateFn_LayoffCostsDivWageBill = @(aprime_val,a_val,z_val,AgentDistMass, tau, w) (-tau*(aprime_val-a_val)*(aprime_val<a_val))/(w*aprime_val); % w*aprime_val is the wage bill
FnsToEvaluateParamNames(5).Names={'p','alpha','cf'};
FnsToEvaluateFn_Production = @(aprime_val,a_val,z_val,AgentDistMass, p, alpha,cf) p*z_val*(aprime_val^alpha)-p*cf;
FnsToEvaluateParamNames(6).Names={'p','alpha'};
FnsToEvaluateFn_Productivity = @(aprime_val,a_val,z_val,AgentDistMass,p,alpha) (p*z_val*(aprime_val^alpha))/(p*aprime_val); % Labor productivity (have delibrately left price p on top and bottom)
FnsToEvaluate={FnsToEvaluateFn_Emp, FnsToEvaluateFn_Hiring, FnsToEvaluateFn_Firing, FnsToEvaluateFn_LayoffCostsDivWageBill, FnsToEvaluateFn_Production, FnsToEvaluateFn_Productivity};

if ImposeFootnote5==1 % Need to deal with the 'special value' in a_grid for new entrants
    % Is simply a matter of replacing a_val with a_val*(a_val~=10^6)
    FnsToEvaluateFn_Hiring = @(aprime_val,a_val,z_val,AgentDistMass) (aprime_val-a_val*(a_val~=10^6))*(aprime_val>a_val*(a_val~=10^6)); % Hiring (need to add the 'firm entry' which involves hiring a single worker)
    FnsToEvaluateFn_Firing = @(aprime_val,a_val,z_val,AgentDistMass) -(aprime_val-a_val*(a_val~=10^6))*(aprime_val<a_val*(a_val~=10^6)); % Firing (need to add the 'firm exits' which involve firing all remaing workers)
    FnsToEvaluateFn_LayoffCostsDivWageBill = @(aprime_val,a_val,z_val,AgentDistMass, tau, w) (-tau*(aprime_val-a_val*(a_val~=10^6))*(aprime_val<a_val*(a_val~=10^6)))/(w*aprime_val); % w*aprime_val is the wage bill
    FnsToEvaluateFn_Production = @(aprime_val,a_val,z_val,AgentDistMass, p, alpha,cf) p*z_val*(aprime_val^alpha)-p*cf*(a_val~=10^6); % the *(a_val~=10^6) relates to footnote 5 on pg 922 that new entrants don't pay the fixed cost of production
    FnsToEvaluate={FnsToEvaluateFn_Emp, FnsToEvaluateFn_Hiring, FnsToEvaluateFn_Firing, FnsToEvaluateFn_LayoffCostsDivWageBill, FnsToEvaluateFn_Production, FnsToEvaluateFn_Productivity};
end

% We will want the aggregate values of these. 
AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
% Entry costs (aggregate)
AggEntryCosts=Params.Ne*Params.p*Params.ce;

output.consumption=(AggValues(5)-AggEntryCosts)/Params.p; % Divide by p to convert nominal to real consumption.

output.averageproductivity=((AggValues(5)-AggEntryCosts)/Params.p)/AggValues(1); % Am guessing this is labor productivity based on (average output)/(average employment) [This calculation includes the fixed-cost of entry] [Division of numerator and denominator by mass cancel out.]
% output.averageproductivity=AggValues(6)/StationaryDist.mass; % Am guessing this is labor productivity at firm level, then averaged (rather than (average output)/(average employment)

output.totalemployment=AggValues(1);

% Since the equilibrium is stationary, and the household is representative,
% the 'utility-adjusted consumption' (more commonly called the 'consumption-equivalent variation') 
% can be calculated from the utility, so just output that.
output.utility_period=log(output.consumption)-Params.a*AggValues(1); % log(c)-aN is the period utility fn of representative household.

% Average Firm Size (i.e., Average number of employees)
output.AvgFirmSize=AggValues(1)/StationaryDist.mass;

output.layoffcostsdivwagebill=AggValues(4)/StationaryDist.mass;
% In stationary eqm, firing must equal hiring, so can use either for
% turnover. (might need to adjust for entry???)
output.TurnoverRateOfJobs=AggValues(2)/AggValues(1); % the "/StationaryDist.mass" cancels top and bottom

% We need a simulated panel based on whole distributions (for calculating
% variance of growth rates and serial correlation in log(n); for survivors).
% Note that because of these two moments we want to calculate it makes more
% sense to have a very large number of two period simulations, and since we
% just want survivors, we won't want entrants.
FnsToEvaluate={FnsToEvaluateFn_Emp, FnsToEvaluateFn_Hiring, FnsToEvaluateFn_Firing};
simoptions.entryinpanel=0; % Don't want entry in this panel data simulation (we are just interested in 'survivors')
simoptions.simperiods=2;
simoptions.numbersims=10^4;
SimPanel=SimPanelValues_Case1(StationaryDist,Policy,FnsToEvaluate,FnsToEvaluateParamNames,Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions, EntryExitParamNames);
Survive_indicator=~isnan(shiftdim((SimPanel(2,2,:)),1));
SimPanel_Survivors=SimPanel(:,:,Survive_indicator);
GrowthRateEmploy=(SimPanel_Survivors(1,2,:)-SimPanel_Survivors(1,1,:))./SimPanel_Survivors(1,1,:);
output.VarianceOfGrowthRate_survivors=var(GrowthRateEmploy);
output.SerialCorrelationLogn_survivors=corr(log(shiftdim(SimPanel_Survivors(1,2,:),2)),log(shiftdim(SimPanel_Survivors(1,1,:),2)));

%% Table 4 stats

% Calculate n_l, the lowest value for which a_prime remains equal to a_val,
% and n_u, the highest value for which a_prime remains equal to a_val.
% That is, for which employment remains equal to last period employment.

% Note that this calculation only makes sense when tau is greater than
% zero. Otherwise there are no adjustment costs of employment, and so no
% 'range' in which we should expect employment to remain unchanged.

if Params.tau>0
    % First, calculate a_prime-a_val: since just interested in when they are equal can actually
    % just do this by first finding the grid point at which the "policy grid index" is equal to current grid index,
    % then report the value of the grid at that index.
    EmploymentDecision_index=shiftdim(Policy,1);
    if ImposeFootnote5==1
        % Need to eliminate the employment decisions of the new entrants
        % as Footnote 5 makes these fundamentaly different.
        EmploymentDecision_index=shiftdim(Policy(1,1:end-1,:),1);
    end
    % Get the first and last 'zeros', for each productivity level
    n_l=nan(n_z,1);
    n_u=nan(n_z,1);
    for z_c=1:n_z
        if max(EmploymentDecision_index(:,z_c))==0
            % Just leave 'nan' if everyone just exits
        else
%             if ImposeFootnote5==0
                n_l(z_c)=a_grid(find(EmploymentDecision_index(:,z_c)==EmploymentDecision_index(1,z_c),1,'last'));
                n_u(z_c)=a_grid(find(EmploymentDecision_index(:,z_c)==EmploymentDecision_index(end,z_c),1,'first'));
%             else
%                 a_grid_temp=a_grid; a_grid_temp(end)=0;
%                 n_l(z_c)=a_grid_temp(find(EmploymentDecision_index(:,z_c)==EmploymentDecision_index(1,z_c),1,'last'));
%                 n_u(z_c)=a_grid_temp(find(EmploymentDecision_index(:,z_c)==EmploymentDecision_index(end,z_c),1,'first'));
%             end
        end
    end    
    output.n_l=n_l;
    output.n_u=n_u;
end

%% Table 5 stats

FnsToEvaluateParamNames(1).Names={'p','alpha'};
FnsToEvaluateFn_MPL = @(aprime_val,a_val,z_val,AgentDistMass,p,alpha) alpha*p*z_val*(aprime_val^(alpha-1)); % MPL
FnsToEvaluate={FnsToEvaluateFn_MPL};
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

NominalMPLvalues=shiftdim(ValuesOnGrid,1);
% Now calcuate the absolute deviations from MPL=1/p as a percentage.
AbsDevsAsPercentage=100*(abs(NominalMPLvalues/Params.p-1/Params.p)/(1/Params.p));

% Partition (a,s) by size of the absolute deviations from MPL=1/p (as %)
FirstPartition=logical((AbsDevsAsPercentage<3));
SecondPartition=logical((AbsDevsAsPercentage>=3).*(AbsDevsAsPercentage<5));
ThirdPartition=logical((AbsDevsAsPercentage>=5).*(AbsDevsAsPercentage<10));
FourthPartition=logical((AbsDevsAsPercentage>=10).*(AbsDevsAsPercentage<15));
FifthPartition=logical((AbsDevsAsPercentage>=15));

% Fraction of firm in each partition (conditional on not exiting)
FractionOfFirmsPerPartition=zeros(5,1);
pdfoffirms=StationaryDist.pdf;
% (Normalized) Mass of Exits
FractionOfExit=sum(sum(pdfoffirms(logical(ExitPolicy))));
FractionOfFirmsPerPartition(1)=sum(sum(pdfoffirms(FirstPartition)));
FractionOfFirmsPerPartition(2)=sum(sum(pdfoffirms(SecondPartition)));
FractionOfFirmsPerPartition(3)=sum(sum(pdfoffirms(ThirdPartition)));
FractionOfFirmsPerPartition(4)=sum(sum(pdfoffirms(FourthPartition)));
FractionOfFirmsPerPartition(5)=sum(sum(pdfoffirms(FifthPartition)));
% Scale for 'conditional on not exiting
FractionOfFirmsPerPartition=FractionOfFirmsPerPartition/(1-FractionOfExit);

output.AbsDevsForTable5=FractionOfFirmsPerPartition;

% Not needed, but for interest/debugging also keep
output.AbsDevsAsPercentage=AbsDevsAsPercentage;
output.NominalMPLdist=NominalMPLvalues;


end
