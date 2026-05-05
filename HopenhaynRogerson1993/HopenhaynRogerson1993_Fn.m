function output=HopenhaynRogerson1993_Fn(ImposeFootnote5, Params,n_a, n_z,a_grid,z_grid,pi_z, ReturnFn, ReturnFnParamNames, DiscountFactorParamNames, EntryExitParamNames, vfoptions, simoptions, heteroagentoptions)
% Largely just a copy-paste of main HopenhaynRogerson1993 code, except
% towards the end when it contains all the calculations of model output
% required for Tables 3, 4 and 5.
% Only substantial difference is that here we treat 'ce' (fixed cost of
% entry) as a parameter, and 'p' as a general equilibrium price to be
% calculated.

n_d=0; % None.
d_grid=[];

% Note: With entry-exit the mass of the distribution of agents often matters. 
% It is possible to use the 'agentmass' extra input argument in all functions to be evaluated.
% If you want to use 'agentmass' as an input to FnsToEvaluate, it must be after the action space but before any of the parameters
FnsToEvaluate.Y = @(aprime,a,z,agentmass,alpha) z*(aprime^alpha); % Real output
% Note: agentmass does nothing in Y, just put it there to demonstate how it can be used.

GEPriceParamNames={'p','Ne'};
GeneralEqmEqns.determineprice = @(RealOutput,p,A) A/RealOutput-p; % The requirement that the price is determined by the demand eqn (or equivalently, can think of this as goods market clearance). You can derive it from FOCs of standard consumption-leisure problem [it is the -U_c/U_N=p/w condition you often see in household problems; remember normalize w=1]: max_{c,N} log(c)-AN s.t. pC=wN+T
GeneralEqmEqns.FreeEntry = @(EValueFn,ce,p,beta) beta*EValueFn-p*ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
if ImposeFootnote5==1
    GeneralEqmEqns.FreeEntry = @(EValueFn,ce,p) EValueFn-p*ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.
end

n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GECondn]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);

Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;

[V,Policy,ExitPolicy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.zeta=1-ExitPolicy;
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);

% The following four are not required for anything, but will output them anyway
output.V=V;
output.Policy=Policy;
output.ExitPolicy=ExitPolicy;
output.StationaryDist=StationaryDist;

%% Now that the stationary equilibrium has been found, generate the statistics required for Tables 3, 4 and 5.

%% Table 3 stats
output.price=Params.p;

FnsToEvaluate.Employment = @(aprime,a,z) aprime; % Employment
FnsToEvaluate.Hiring = @(aprime,a,z) (aprime-a)*(aprime>a); % Hiring (need to add the 'firm entry' which involves hiring a single worker)
FnsToEvaluate.Firing = @(aprime,a,z) -(aprime-a)*(aprime<a); % Firing (need to add the 'firm exits' which involve firing all remaing workers)
FnsToEvaluate.LayoffCostsDivWageBill = @(aprime,a,z, tau, w) (-tau*(aprime-a)*(aprime<a))/(w*aprime); % w*aprime is the wage bill
FnsToEvaluate.Production = @(aprime,a,z, p, alpha,cf) p*z*(aprime^alpha)-p*cf;
FnsToEvaluate.Productivity = @(aprime,a,z,p,alpha) (p*z*(aprime^alpha))/(p*aprime); % Labor productivity (have delibrately left price p on top and bottom)

if ImposeFootnote5==1 % Need to deal with the 'special value' in a_grid for new entrants
    % Is simply a matter of replacing a with a*(a~=10^6)
    FnsToEvaluate.Hiring = @(aprime,a,z) (aprime-a*(a~=10^6))*(aprime>a*(a~=10^6)); % Hiring (need to add the 'firm entry' which involves hiring a single worker)
    FnsToEvaluate.Firing = @(aprime,a,z) -(aprime-a*(a~=10^6))*(aprime<a*(a~=10^6)); % Firing (need to add the 'firm exits' which involve firing all remaing workers)
    FnsToEvaluate.LayoffCostsDivWageBill = @(aprime,a,z, tau, w) (-tau*(aprime-a*(a~=10^6))*(aprime<a*(a~=10^6)))/(w*aprime); % w*aprime is the wage bill
    FnsToEvaluate.Production = @(aprime,a,z, p, alpha,cf) p*z*(aprime^alpha)-p*cf*(a~=10^6); % the *(a~=10^6) relates to footnote 5 on pg 922 that new entrants don't pay the fixed cost of production
end

% We will want the aggregate values of these. 
AggVars=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions,EntryExitParamNames);
% Entry costs (aggregate)
AggEntryCosts=Params.Ne*Params.p*Params.ce;

output.consumption=(AggVars.Production.Mean-AggEntryCosts)/Params.p; % Divide by p to convert nominal to real consumption.

output.averageproductivity=((AggVars.Production.Mean-AggEntryCosts)/Params.p)/AggVars.Employment.Mean; % Am guessing this is labor productivity based on (average output)/(average employment) [This calculation includes the fixed-cost of entry] [Division of numerator and denominator by mass cancel out.]
% output.averageproductivity=AggValues.Productivity.Mean/StationaryDist.mass; % Am guessing this is labor productivity at firm level, then averaged (rather than (average output)/(average employment)

output.totalemployment=AggVars.Employment.Mean;

% Since the equilibrium is stationary, and the household is representative,
% the 'utility-adjusted consumption' (more commonly called the 'consumption-equivalent variation') 
% can be calculated from the utility, so just output that.
output.utility_period=log(output.consumption)-Params.a*AggVars.Employment.Mean; % log(c)-aN is the period utility fn of representative household.

% Average Firm Size (i.e., Average number of employees)
output.AvgFirmSize=AggVars.Employment.Mean/StationaryDist.mass;

output.layoffcostsdivwagebill=AggVars.LayoffCostsDivWageBill.Mean/StationaryDist.mass;
% In stationary eqm, firing must equal hiring, so can use either for
% turnover. (might need to adjust for entry???)
output.TurnoverRateOfJobs=AggVars.Hiring.Mean/AggVars.Employment.Mean; % the "/StationaryDist.mass" cancels top and bottom

% We need a simulated panel based on whole distributions (for calculating
% variance of growth rates and serial correlation in log(n); for survivors).
% Note that because of these two moments we want to calculate it makes more
% sense to have a very large number of two period simulations, and since we
% just want survivors, we won't want entrants.
FnsToEvaluate2.Employment=FnsToEvaluate.Employment;
FnsToEvaluate2.Hiring=FnsToEvaluate.Hiring;
FnsToEvaluate2.Firing=FnsToEvaluate.Firing;
simoptions.entryinpanel=0; % Don't want entry in this panel data simulation (we are just interested in 'survivors')
simoptions.simperiods=2;
simoptions.numbersims=10^4;
SimPanel=SimPanelValues_InfHorz(StationaryDist,Policy,FnsToEvaluate2,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z, simoptions, EntryExitParamNames);
Survive_indicator=~isnan(SimPanel.Employment(2,:));
SimPanel_Survivors.Employment=SimPanel.Employment(:,Survive_indicator);
SimPanel_Survivors.Hiring=SimPanel.Hiring(:,Survive_indicator);
SimPanel_Survivors.Firing=SimPanel.Firing(:,Survive_indicator);
GrowthRateEmploy=(SimPanel_Survivors.Employment(2,:)-SimPanel_Survivors.Employment(1,:))./SimPanel_Survivors.Employment(1,:);
output.VarianceOfGrowthRate_survivors=var(GrowthRateEmploy);
output.SerialCorrelationLogn_survivors=corr(log(SimPanel_Survivors.Employment(2,:)'),log(SimPanel_Survivors.Employment(1,:)'));

%% Table 4 stats

% Calculate n_l, the lowest value for which a_prime remains equal to a,
% and n_u, the highest value for which a_prime remains equal to a.
% That is, for which employment remains equal to last period employment.

% Note that this calculation only makes sense when tau is greater than
% zero. Otherwise there are no adjustment costs of employment, and so no
% 'range' in which we should expect employment to remain unchanged.

if Params.tau>0
    % First, calculate a_prime-a: since just interested in when they are equal can actually
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

FnsToEvaluate3.MPL = @(aprime,a,z,p,alpha) alpha*p*z*(aprime^(alpha-1)); % MPL
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy, FnsToEvaluate3, Params, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions,EntryExitParamNames,StationaryDist); % StationaryDist is only needed because using entry-exit

% Now calcuate the absolute deviations from MPL=1/p as a percentage.
AbsDevsAsPercentage=100*(abs(ValuesOnGrid.MPL/Params.p-1/Params.p)/(1/Params.p));

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
output.NominalMPLdist=ValuesOnGrid.MPL;


end
