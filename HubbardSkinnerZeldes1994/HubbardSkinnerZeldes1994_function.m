function [Table1row, Table2, Table3, AgeConditionalStats]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,Names_i,simoptions, PTypeDistParamNames)

%% This code replicates the results of 
% Hubbard, Skinner & Zeldes (1994) - The importance of precautionary motives in explaining individual and aggregate saving
% The model detailed in this Carnegie-Rochester Conference Series on Public Policy paper appears to also be exactly the same as that used in
% Hubbard, Skinner & Zeldes (1994) - Expanding the Life-Cycle Model: Precautionary Saving and Public Policy
% Hubbard, Skinner & Zeldes (1994) - Precautionary Saving and Social Insurance
% published in AER:P&P and JPE respecitively. Although I do not replicate all the results of these other two.

% Note that the model includes a floor on consumption. This creates non-convexities in the
% intertemporal budget constraint and means that basic Euler eqn based methods would not work.
% HSZ use an Euler eqn based method, but have to add considerable complications in order for it 
% to still work (see pg 113). This makes it a nice example for our present purposes as discertized VFI methods 
% of the VFI Toolkit are able to handle such models with no extra effort required.

%% State space variables of the model
% Age
% One endogenous variable: assets
% Two stochastic exogenous variables: income shock, medical expense shock
% Three permanent types: 
%  pt1: high-school dropout
%  pt2: high-school
%  pt3: university

%%

% Params.J=80; % Ages 21 to 100 inclusive.
% WorkingAgeVec=21:1:65;
% RetiredAgeVec=66:100;
% Params.age=[WorkingAgeVec,RetiredAgeVec];

% 
% n_a=250;
maxa=5*10^5;
% n_z=[15,15]; % income, medical
% Names_i={'NoHighSchool','HighSchool','College'}; % Number of fixed types
N_j=Params.J; % Number of periods in finite horizon
% Params.q=3; % For tauchen method. They used q=2.5: pg 112 of pdf, "range between 2.5 standard deviations (of the unconditional distribution) above and below..."
% 
% Params.gamma=3;
% % Gamma plays three roles:      
%         % gamma is coefficient of relative risk aversion
%         % 1/gamma is the intertemporal elasticity of substitution of consumption
%         % gamma+1 is the coefficient of relative prudence (Kimball, 1990)
% 
% Params.delta=0.03; % rate at which agents discount future
Params.beta=1/(1+Params.delta);
% Params.r=0.03; % interest rate on assets
% 
% % Mortality rates (probability of survival)
% % Table A.1
% Params.dj=[0.00058, 0.00061, 0.00062, 0.00064, 0.00065, 0.00067, 0.00069, 0.00070, 0.00072, 0.00075, 0.00078, 0.00082, 0.00086, 0.00091, 0.00098, 0.00105, 0.00115, 0.00128, 0.00144, 0.00161, 0.00180, 0.00200,...
%     0.00221, 0.00242, 0.00266, 0.00292, 0.00320, 0.00349, 0.00380, 0.00413, 0.00450, 0.00490, 0.00533,...
%     0.00533, 0.00581, 0.00632, 0.00689, 0.00749, 0.00811, 0.00878, 0.00952, 0.01033, 0.01121, 0.01223,...
%     0.01332, 0.01455, 0.01590, 0.01730, 0.01874, 0.02028, 0.02203, 0.02404, 0.02623, 0.02863, 0.03128,...
%     0.03432, 0.03778, 0.04166, 0.04597, 0.05078, 0.05615, 0.06214, 0.06885, 0.07631, 0.08455, 0.09352,...
%     0.10323, 0.11367, 0.12484, 0.13677, 0.14938, 0.16289, 0.17721, 0.19234, 0.20828, 0.22418, 0.23980, 0.25495, 0.26937, 0.28284]; % conditional probability of death
% Params.sj=1-Params.dj; % Conditional survival probabilities. Act as a discount rate.
% 
% % Wage (Earnings) incomes, W
% % Table A.2
% Params.DeterministicWj_working.NoHighSchool=[10993-196*WorkingAgeVec+2167*((WorkingAgeVec.^2)/100)-2655*((WorkingAgeVec.^3)/10000)];
% Params.DeterministicWj_working.HighSchool=[-11833+1421*WorkingAgeVec-137*((WorkingAgeVec.^2)/100)-2186*((WorkingAgeVec.^3)/10000)]; 
% Params.DeterministicWj_working.College=[72270-5579*WorkingAgeVec+19200*((WorkingAgeVec.^2)/100)-18076*((WorkingAgeVec.^3)/10000)];
% % Table A.3
% Params.DeterministicWj_retired.NoHighSchool=[17861-133*RetiredAgeVec]; 
% Params.DeterministicWj_retired.HighSchool=[29733-245*RetiredAgeVec]; 
% Params.DeterministicWj_retired.College=[48123-429*RetiredAgeVec];
% Params.DeterministicWj.NoHighSchool=[Params.DeterministicWj_working.NoHighSchool, Params.DeterministicWj_retired.NoHighSchool];
% Params.DeterministicWj.HighSchool=[Params.DeterministicWj_working.HighSchool, Params.DeterministicWj_retired.HighSchool];
% Params.DeterministicWj.College=[Params.DeterministicWj_working.College, Params.DeterministicWj_retired.College];
% % Compare to Fig A.1(a-c): they won't be exactly the same as drop the year
% % fixed effects, but eyeballing suggests they are fine.
% plot(21:1:100, Params.DeterministicWj.NoHighSchool, 21:1:100, Params.DeterministicWj.HighSchool, 21:1:100, Params.DeterministicWj.College)
% % Stochastic Wj (Table A.4): u_it, an AR(1) in regressions in paper
% Params.w_rho=[0.955; 0.946; 0.955];
% Params.w_sigmasqepsilon=[0.033; 0.025; 0.016];
% Params.w_sigmasqu=Params.w_sigmasqepsilon./(1-Params.w_rho.^2);
% Params.w_sigmasqupsilon=[0.040; 0.021; 0.014]; % Estimated from PSID but not used in model.
[z1_grid.NoHighSchool,pi_z1.NoHighSchool]=discretizeAR1_Tauchen(0,Params.w_rho(1),sqrt(Params.w_sigmasqepsilon(1)),n_z(1),Params.q);
[z1_grid.HighSchool,pi_z1.HighSchool]=discretizeAR1_Tauchen(0,Params.w_rho(2),sqrt(Params.w_sigmasqepsilon(2)),n_z(1),Params.q);
[z1_grid.College,pi_z1.College]=discretizeAR1_Tauchen(0,Params.w_rho(3),sqrt(Params.w_sigmasqepsilon(3)),n_z(1),Params.q);
% % Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% % "This specification ensures that when we compare the certainty case with
% % the earnings uncertainty case, we hold the age-conditional means of earnings constant."
% 
% % Medical Expenses, M
% % Table A.5
% Params.DeterministicMj_working=[5.373+0.073*WorkingAgeVec-0.753*((WorkingAgeVec.^2)/1000); 4.749+0.106*WorkingAgeVec-1.084*((WorkingAgeVec.^2)/1000); 4.816+0.109*WorkingAgeVec-1.090*((WorkingAgeVec.^2)/1000)];
% Params.DeterministicMj_retired=[14.441-0.200*RetiredAgeVec+1.288*((RetiredAgeVec.^2)/1000); 11.371-0.101*RetiredAgeVec+0.540*((RetiredAgeVec.^2)/1000); 9.553-0.054*RetiredAgeVec-0.297*((RetiredAgeVec.^2)/1000)];
% Params.DeterministicMj=[Params.DeterministicMj_working, Params.DeterministicMj_retired]; % is in logs
% Params.m_rho=0.901*ones(3,1);
% Params.m_sigmasqepsilon=[0.175; 0.156; 0.153];
% Params.m_sigmasqomega=[0.220; 0.220; 0.220];
[z2_grid.NoHighSchool,pi_z2.NoHighSchool]=discretizeAR1_Tauchen(0,Params.m_rho(1),sqrt(Params.m_sigmasqepsilon(1)),n_z(2),Params.q);
[z2_grid.HighSchool,pi_z2.HighSchool]=discretizeAR1_Tauchen(0,Params.m_rho(2),sqrt(Params.m_sigmasqepsilon(2)),n_z(2),Params.q);
[z2_grid.College,pi_z2.College]=discretizeAR1_Tauchen(0,Params.m_rho(3),sqrt(Params.m_sigmasqepsilon(3)),n_z(2),Params.q);

% % Consumption Floor
% Params.Cbar=7000; % (middle of pg. 111)

%% Grids
z_grid.NoHighSchool=[z1_grid.NoHighSchool; z2_grid.NoHighSchool];
z_grid.HighSchool=[z1_grid.HighSchool; z2_grid.HighSchool];
z_grid.College=[z1_grid.College; z2_grid.College];

pi_z.NoHighSchool=kron(pi_z2.NoHighSchool,pi_z1.NoHighSchool); % note, kron() in reverse order
pi_z.HighSchool=kron(pi_z2.HighSchool,pi_z1.HighSchool); % note, kron() in reverse order
pi_z.College=kron(pi_z2.College,pi_z1.College); % note, kron() in reverse order

a_grid=linspace(0,maxa,n_a)'; % Could probably do better by adding more grid near zero

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(aprime,a,W_z1,M_z2,age,gamma,r,Cbar,DeterministicWj, w_sigmasqu, DeterministicMj, m_sigmasqmew) HubbardSkinnerZeldes1994_ReturnFn(aprime,a,W_z1,M_z2,age,gamma,r,Cbar,DeterministicWj, w_sigmasqu, DeterministicMj, m_sigmasqmew)

%% Now solve the value function iteration problem

vfoptions.verbose=1;
tic;
% Don't keep V as it is not needed for anything
[~, Policy]=ValueFnIter_Case1_FHorz_PType(0,n_a,n_z,N_j,Names_i, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,vfoptions);
toc

% %% Draw some policy functions of 'fixed type 1'.
% % These are not part of replications, just to illustrate some uses of VFI Toolkit.
% % Note that what is being plotted are 'policy indexes', not 'policy values'.
% figure(2)
% surf(reshape(Policy.NoHighSchool(1,:,8,8,:),[250,80])-kron(linspace(1,250,250)',ones(1,80)))
% figure(3)
% surf(reshape(Policy.NoHighSchool(1,:,15,8,:),[250,80])-kron(linspace(1,250,250)',ones(1,80)))
% figure(4)
% surf(reshape(Policy.NoHighSchool(1,:,1,8,:),[250,80])-kron(linspace(1,250,250)',ones(1,80)))

%% Create some outputs to replicate those of HubbardSkinnerZeldes1994
% 
% % Note: The annual growth rate of the population is assumed to be 1%.
% Params.n=0.01;
% 
% simoptions.numbersims=10^3; % Number of simulations on which panel data (and life-cycle profile) results will be based.

%% Create an initial distribution from which agents are drawn/born.

% I have been unable to find any mention in the paper of how the initial
% distribution from which agents are born is determined. I therefore simply
% assume that they are all born with zero assets, the stationary distribution 
% on income shocks, and the mean shock value on medical expenses.
% (Paper does specify that all 'newborn' are age 21; ie. first period.)
% More tricky is what weight to attach to each of the permanent types, this
% remains unclear. HubbardSkinnerZeldes1994 do not appear to report these
% numbers in the paper at all.

z1staty.NoHighSchool=ones(size(z1_grid.NoHighSchool))/length(z1_grid.NoHighSchool);
for ii=1:1000
    z1staty.NoHighSchool=pi_z1.NoHighSchool'*z1staty.NoHighSchool;
end
z1staty.NoHighSchool=z1staty.NoHighSchool./sum(z1staty.NoHighSchool(:)); % Normalize to 1 (I was getting numerical rounding errors)
z1staty.HighSchool=ones(size(z1_grid.HighSchool))/length(z1_grid.HighSchool);
for ii=1:1000
    z1staty.HighSchool=pi_z1.HighSchool'*z1staty.HighSchool;
end
z1staty.HighSchool=z1staty.HighSchool./sum(z1staty.HighSchool(:)); % Normalize to 1 (I was getting numerical rounding errors)
z1staty.College=ones(size(z1_grid.College))/length(z1_grid.College);
for ii=1:1000
    z1staty.College=pi_z1.College'*z1staty.College;
end
z1staty.College=z1staty.College./sum(z1staty.College(:)); % Normalize to 1 (I was getting numerical rounding errors)

InitialDist.NoHighSchool=zeros([n_a,n_z]); 
InitialDist.NoHighSchool(1,:,ceil(n_z(2)/2))=z1staty.NoHighSchool';
InitialDist.HighSchool=zeros([n_a,n_z]); 
InitialDist.HighSchool(1,:,ceil(n_z(2)/2))=z1staty.HighSchool';
InitialDist.College=zeros([n_a,n_z]); 
InitialDist.College(1,:,ceil(n_z(2)/2))=z1staty.College';
% Note: We just put each permanent type in as a mass of one. The toolkit
% knows that it should then combine this with the PTypeDistParamNames entry
% in the Params structure for the relative masses of each permanent type.
% To put it another way, InitialDist.College contains the distribution of
% agents conditional on permanent type being College.

% Calculate the relative masses of the different ages
% Taking into account population growth n and (conditional) survival rates sj.
ageweights=cumprod(Params.sj).*cumprod((1/(1+Params.n))*ones(1,Params.J)); % First part is based on survival, second part on population growth
ageweights=ageweights./sum(ageweights); % Normalize to sum to 1
Params.mewj=ageweights;
AgeWeightsParamNames={'mewj'};


%% Life-cycle profiles

FnsToEvaluate.Assets = @(aprime,a,z1,z2) a; 
FnsToEvaluate.Earnings= @(aprime,a,z1,z2,DeterministicWj,w_sigmasqu) exp(log(DeterministicWj)-0.5*w_sigmasqu+z1);
FnsToEvaluate.age = @(aprime,a,z1,z2,age) age;
FnsToEvaluate.Income = @(aprime,a,z1,z2,r,DeterministicWj,w_sigmasqu) r*a+exp(log(DeterministicWj)-0.5*w_sigmasqu+z1);  
FnsToEvaluate.Consumption= @(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj) HubbardSkinnerZeldes1994_ConsumptionFn(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj);
FnsToEvaluate.TR = @(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj) HubbardSkinnerZeldes1994_TRFn(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj);
FnsToEvaluate.AssetIncomeRatio = @(aprime,a,z1,z2,r,DeterministicWj,w_sigmasqu) a/(r*a+exp(log(DeterministicWj)-0.5*w_sigmasqu+z1));  
FnsToEvaluate.SavingsRate = @(aprime,a,z1,z2,r,DeterministicWj,w_sigmasqu) (aprime-a)/(r*a+exp(log(DeterministicWj)-0.5*w_sigmasqu+z1));  
FnsToEvaluate.GrossSavings = @(aprime,a,z1,z2) aprime-a;  
FnsToEvaluate.ptype = @(aprime,a,z1,z2,hhtype) hhtype;

fprintf('Doing StationaryDist and AgeConditionalStats \n')
% the command for the life-cycle profiles needs the agent stationary distribution, so we do that first
StationaryDist=StationaryDist_Case1_FHorz_PType(InitialDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,0,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy, FnsToEvaluate, Params,0,n_a,n_z,N_j,Names_i,[],a_grid,z_grid, simoptions);

% % Figure: Assets
% figure(1)
% plot(Params.age,AgeConditionalStats.Assets.NoHighSchool.Mean,Params.age,AgeConditionalStats.Assets.HighSchool.Mean,Params.age,AgeConditionalStats.Assets.College.Mean)
% % Figure: Earnings
% figure(2)
% plot(Params.age,AgeConditionalStats.Assets.NoHighSchool.Mean,Params.age,AgeConditionalStats.Assets.HighSchool.Mean,Params.age,AgeConditionalStats.Assets.College.Mean)

% Simulate Panel Data
fprintf('Simulate Panel Data \n')
% Same variables as we used for the life-cycle profiles.
SimPanelValues=SimPanelValues_FHorz_Case1_PType(InitialDist,PTypeDistParamNames,Policy, FnsToEvaluate,Params,0,n_a,n_z,N_j,Names_i,0,a_grid,z_grid,pi_z, simoptions);

fprintf('Put results together for this parameter combination \n')


%% Table 1: Asset-Income ratio and Savings Rate (aggregate and conditional on fixed-type)

% At first, not clear if Table 1 reports the mean of the ratio, or the ratio of the mean (analagously for savings rate)
% From Appendix B it is clear that Table 1 reports the ratio of the mean.
% Following is based on Appendix B: Constructing aggregate consumption, earnings, and assets. In the sense of following 
% the variable definitions there

% This is pretty trivial, we just use AggVars
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,0,n_a,n_z,N_j,Names_i,[],a_grid,z_grid,simoptions);

AssetIncomeRatioByFixedType=[AggVars.Assets.NoHighSchool.Mean/AggVars.Income.NoHighSchool.Mean,...
    AggVars.Assets.HighSchool.Mean/AggVars.Income.HighSchool.Mean,...
    AggVars.Assets.College.Mean/AggVars.Income.College.Mean];
AssetIncomeRatio_Agg=AggVars.Assets.Mean/AggVars.Income.Mean;

SavingsRateByFixedType=[AggVars.GrossSavings.NoHighSchool.Mean/AggVars.Income.NoHighSchool.Mean,...
    AggVars.GrossSavings.HighSchool.Mean/AggVars.Income.HighSchool.Mean,...
    AggVars.GrossSavings.College.Mean/AggVars.Income.College.Mean];
SavingsRate_Agg=AggVars.GrossSavings.Mean/AggVars.Income.Mean;

% Income here is after tax but before transfers.

Table1row=[Params.delta, Params.gamma, AssetIncomeRatioByFixedType, AssetIncomeRatio_Agg, SavingsRateByFixedType, SavingsRate_Agg];
% Note: Hubbard, Skinner & Zeldes (1994) calculate all of these from the panel data, but VFI Toolkit has commands that calculate these 
% more directly so we have used those. Calculating them from panel data is not difficult, but just requires writing a whole bunch of 
% unnecessary code.


%% Table 2: creates the numbers relevant for Table 2
% Conditional probabilities of having consumption approximately equal to income.

HighCorrIncomeCons=zeros(6,3);
NumberOfHHs=zeros(6,3);
LowAssets_NumberOfHHs=zeros(6,3);
LowAssets_HighCorrIncomeCons=zeros(6,3);

% AverageIncomeByAgeBin=[mean(AgeConditionalStats(4,1:9,1)),mean(AgeConditionalStats(4,10:19,1)),mean(AgeConditionalStats(4,20:29,1)),mean(AgeConditionalStats(4,30:39,1)),mean(AgeConditionalStats(4,40:49,1)),mean(AgeConditionalStats(4,50:end,1))];

for ii=1:simoptions.numbersims
    for jj=1:80
        cons=SimPanelValues.Consumption(jj,ii);
        income=SimPanelValues.Income(jj,ii);
        assets=SimPanelValues.Assets(jj,ii);
        age=SimPanelValues.age(jj,ii);
        hhtype=SimPanelValues.ptype(jj,ii);
        
        if age<=29
            agebin=1;
        elseif age<=39
            agebin=2;
        elseif age<=49
            agebin=3;
        elseif age<=59
            agebin=4;
        elseif age<=69
            agebin=5;
        else
            agebin=0;
        end
        
        % Do the needed calculations for each agebin
        if agebin>0
            NumberOfHHs(agebin,hhtype)=NumberOfHHs(agebin,hhtype)+1;
            if (cons/income)>=0.95 && (cons/income)<=1.05
                HighCorrIncomeCons(agebin,hhtype)=HighCorrIncomeCons(agebin,hhtype)+1;
            end
            
            if assets<income % Is not quite clear from paper which HHs in simulated data are considered to satisfy the 'Initial Assets <0.5*Average Income' criterion.
                LowAssets_NumberOfHHs(agebin,hhtype)=LowAssets_NumberOfHHs(agebin,hhtype)+1;
                if (cons/income)>=0.95 && (cons/income)<=1.05
                    LowAssets_HighCorrIncomeCons(agebin,hhtype)=LowAssets_HighCorrIncomeCons(agebin,hhtype)+1;
                end
            end
        end
        % Also do the Total in the sixth row
        NumberOfHHs(6,hhtype)=NumberOfHHs(6,hhtype)+1;
        if (cons/income)>=0.95 && (cons/income)<=1.05
            HighCorrIncomeCons(6,hhtype)=HighCorrIncomeCons(6,hhtype)+1;
        end
        if assets<income % Is not quite clear from paper which HHs in simulated data are considered to satisfy the 'Initial Assets <0.5*Average Income' criterion.
            LowAssets_NumberOfHHs(6,hhtype)=LowAssets_NumberOfHHs(6,hhtype)+1;
            if (cons/income)>=0.95 && (cons/income)<=1.05
                LowAssets_HighCorrIncomeCons(6,hhtype)=LowAssets_HighCorrIncomeCons(6,hhtype)+1;
            end
        end
    end
end

Table2=[HighCorrIncomeCons./NumberOfHHs;LowAssets_HighCorrIncomeCons./LowAssets_NumberOfHHs];

%% Table 3 
% Campbell-Mankiw-Lusardi Euler Eqns
% My impression from paper is that these regressions from (pooled) simulated panel
% do not correct for the 1% population growth that they assume (ie. that it
% is ignored here). I therefore follow this.
% According to Hubbard, Skinner & Zeldes (1994) the original Cambell-Mankiw
% regressions were on aggregate consumption and income, but in this model
% both of those are constants so this regression would not be possible.
% Presumably it has therefore been done (simulated) microdata instead.

% cons=SimPanelValues(5,:,:);
% income=SimPanelValues(4,:,:);

vecsize=[(N_j-1)*simoptions.numbersims,1];

DeltaC=reshape(SimPanelValues.Consumption(2:end,:)-SimPanelValues.Consumption(1:(end-1),:),vecsize);
DeltaY=reshape(SimPanelValues.Income(2:end,:)-SimPanelValues.Income(1:(end-1),:),vecsize);
DeltalnC=reshape(log(SimPanelValues.Consumption(2:end,:))-log(SimPanelValues.Consumption(1:(end-1),:)),vecsize);
DeltalnY=reshape(log(SimPanelValues.Income(2:end,:))-log(SimPanelValues.Income(1:(end-1),:)),vecsize);
age=reshape(SimPanelValues.age(2:end,:),vecsize);
constant=ones(vecsize);

% All the regressions appear to be two-stage least squared, using lags as
% instruments for the change in income terms.

% First regression
y2=DeltaY(4:end);
X2=[constant(4:end),DeltaC(2:end-2),DeltaC(1:end-3),DeltaY(2:end-2),DeltaY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltaYfitted=X2*b2;
y=DeltaC(4:end);
X=[constant(4:end),DeltaYfitted];
[bcolumn1,bintcolumn1,r,rint,stats]=regress(y,X);
% Second regression
y2=DeltalnY(4:end);
X2=[constant(4:end),DeltalnC(2:end-2),DeltalnC(1:end-3),DeltalnY(2:end-2),DeltalnY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltalnYfitted=X2*b2;
y=DeltalnC(4:end);
X=[constant(4:end),DeltalnYfitted];
[bcolumn2,bintcolumn2,r,rint,stats]=regress(y,X);
% Third regression
y2=DeltalnY(4:end);
X2=[constant(4:end),DeltalnC(3:end-1),DeltalnC(2:end-2),DeltalnC(1:end-3),DeltalnY(3:end-1),DeltalnY(2:end-2),DeltalnY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltalnYfitted=X2*b2;
y=DeltalnC(4:end);
X=[constant(4:end),DeltalnYfitted];
[bcolumn3,bintcolumn3,r,rint,stats]=regress(y,X);
% Fourth regression
y2=DeltalnY(4:end);
X2=[constant(4:end),DeltalnC(3:end-1),DeltalnC(2:end-2),DeltalnC(1:end-3),DeltalnY(3:end-1),DeltalnY(2:end-2),DeltalnY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltalnYfitted=X2*b2;
y=DeltalnC(4:end);
X=[constant(4:end),DeltalnYfitted,age(4:end),age(4:end).^2];
[bcolumn4,bintcolumn4,r,rint,stats]=regress(y,X);

% These regressions ignore the surivival probabilies. I suspect that Hubbard, Skinner & Zeldes (1994) do not ignore these.

Table3=nan(8,4);
Table3(1,1)=bcolumn1(2);
Table3(2,1)=bcolumn1(2)/((bintcolumn1(2,2)-bcolumn1(2))/1.96); % Coefficient estimate divided by standard error. Since matlab just gives the 95% confidence intervals (which correspond to plus and minus 1.96 std errors) we have to calculate the standard error in the coefficient.
Table3(3,2)=bcolumn2(2);
Table3(4,2)=bcolumn2(2)/((bintcolumn2(2,2)-bcolumn2(2))/1.96);
Table3(3,3)=bcolumn3(2);
Table3(4,3)=bcolumn3(2)/((bintcolumn3(2,2)-bcolumn3(2))/1.96);
Table3(3,4)=bcolumn4(2);
Table3(4,4)=bcolumn4(2)/((bintcolumn4(2,2)-bcolumn4(2))/1.96);
Table3(5,4)=bcolumn4(3);
Table3(6,4)=bcolumn4(3)/((bintcolumn4(3,2)-bcolumn4(3))/1.96);
Table3(7,4)=bcolumn4(4);
Table3(8,4)=bcolumn4(4)/((bintcolumn4(4,2)-bcolumn4(4))/1.96);



end
