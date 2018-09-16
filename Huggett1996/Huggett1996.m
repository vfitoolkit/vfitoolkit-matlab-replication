% Replication of Huggett (1996) - Wealth Distribution in Life Cycle Economies
%
% This code takes much longer than necessary to solve as I wish to double-
% and triple-check the step of finding general equilibrium as part of the
% replication.
%
% I follow notation of Huggett (1996) for all parameters except the number
% of period (and retirement age) which I denote by J (JR) rather than N.
%
% One endogenous variable: assets
% One stochastic exogenous variables: income shock
% Age

CheckUniquenessOfGE=0; % If =1, will used grid on prices around the GE to check for possibility of other GE. Major increase in run time.

Params.J=79; % Ages 20 to 98 inclusive.

% Grid sizes to use
n_a=1501;
% n_z=19; % income (these 19 points are hardcoded into z_grid and pi_z, done this way due to how Huggett sets them up)
N_j=Params.J; % Number of periods in finite horizon

% A few lines needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
try % Server has 16 cores, but is shared with other users, so use max of 8.
    parpool(8)
    gpuDevice(1)
catch % Desktop has less than 8, so will give error, on desktop it is fine to use all available cores.
    parpool
end
PoolDetails=gcp;
NCores=PoolDetails.NumWorkers;
simoptions.ncores=NCores;

simoptions.parallel=3;

%% Declare the model parameters
% Note that r, w, and G will be determined in General equilbrium, so these are really just initial guesses.

% Preference
Params.beta=1.011; % rate at which agents discount future (is greater than one due to the 'stochastic probability of death')
Params.sigma=1.5; % Risk-aversion

Params.bc_equalsminusw=0; % An indicator variable for the first of the two possible values for the borrowing constraint, 0 and -w

% Demographics
Params.JR=46; % Retirement at age 65 (Note: never actually used directly as is implicit in the deterministic earnings profile and the retirement benefits
Params.n=0.012; % Population growth rate of 1%

% Tax rates
% Params.tau % Determined based on r: Params.tau=0.195*(1-Params.delta*K/Y) % Note that K/Y can be calculated from r, see below
Params.tau=0.195*(1-0.06*3); % Just give an initial guess for tau here
Params.theta=0.1;
% Accidental bequests
% Params.T % Determined in GE
% Retirement benefits: I set them equal to b*bvec (in Huggetts notation this is just b which is itself a function of retirement status, 
% I seperate it into a scalar and an indicator on whether or not you are retired).
Params.bvec=[zeros(1,Params.JR-1),ones(1,1+Params.J-Params.JR)]; % Set it up as an age dependent vector that is zero before retirement age
% Note: bvec is really just an indicator of retirement

% Production function
Params.A=0.895944;
Params.alpha=0.36; % Capital-share in Cobb-Douglas production function (I currently live in New Zealand where this is the actual capital-share of GDP ;)
Params.delta=0.06; % Depreciation rate

% Survival probabilities, based on Jordan, C., Life contingencies, 2nd ed. (Society of Actuaries).
% Huggett (1996) dates it as the 1975 print, the following link is to the 1991 imprint. But since it is still the 2nd edition (which first appeared 1967)
% Seems to be the right numbers. (the pg 342 of the pdf, 346 of book, appears to be the closest thing)
% https://vdocuments.mx/download/life-contingencies-chester-wallace-jordanpdf
Params.dj=[0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
    0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
    0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
    0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
    0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
    0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
    0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
    0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018]; % conditional probability of death
Params.sj=1-Params.dj; % Conditional survival probabilities. Act as a discount rate.

% Declare the age dependent parameters. This is a simple matter of creating
% the parameter as a row vector of length J (the VFI Toolkit automatically
% detects which parameters depend on age and which do not).
% Params.ybarj=log(0.5289)+log([linspace(0.3,1.25,12),linspace(1.3,1.5,6),linspace(1.5,1.5,11),linspace(1.49,1.05,11),linspace(1,0.1,9),linspace(0.08,0,5),zeros(1,Params.J-54)]); % deterministic income depends on age
Corbae_deterministicincome = [0.0911 0.1573 0.2268 0.2752 0.3218 0.3669 0.4114 0.4559 0.4859 0.5164 0.5474 0.5786 0.6097 0.6311 0.6517 0.6711 0.6893 0.7060 0.7213 0.7355 0.7489 0.7619 0.7747 0.7783 0.7825 0.7874 0.7931 0.7994 0.7923 0.7850 0.7771 0.7679 0.7567 0.7351 0.7105 0.6822 0.6500 0.6138 0.5675 0.5183 0.4672 0.3935 0.3239 0.2596 0.1955 0.1408 0.0959 0.0604 0.0459 0.0342 0.0246 0.0165 0.0091 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % Figure 1 of Huggett (1996) shows: plot(19+(1:1:Params.J),exp(Params.ybarj))
Params.ybarj=log(Corbae_deterministicincome);
% The log(0.5289)+ come from a 2014 Homework Handout by Dean Corbae who
% states that he got them from Mark Huggett (I got the Corbae handout files from a friend).
% The only thing that is odd is the normalization by 0.5289. I can't find
% mention of this in Huggett (1996) but according to the Corbae Homework
% it is about getting the ratio of 'work' L to population N correct; to
% 0.5289 as in data.

% Stochastic income, y:
Params.gamma=0.96;
Params.sigmasqepsilon=0.045;
% Params.sigmasqy=Params.sigmasqepsilon./(1-Params.gamma.^2); 
% Initial distribution of income
Params.sigmasqy1=0.38;

% Huggett use 17 states equally spaced between +-4*sigmasqy1, with an 18th
% state at 6*sigmasqy1, "The transition probabilities between states are
% calculated by integrating the area under the normal distribution conditional on
% the current value of the state."
n_z=18;
% z_grid=[linspace(-4*sqrt(Params.sigmasqy1),4*sqrt(Params.sigmasqy1),17),6*sqrt(Params.sigmasqy1)]';
% pi_z=nan(18,18);
% % Following lines implement the transition matrix, they are largely just a copy of some code from the TauchenMethod() command.
% sigma=sqrt(Params.sigmasqepsilon); %stddev of e
% for ii=1:length(z_grid)
%     pi_z(ii,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2-Params.gamma*z_grid(ii),0,sigma);
%     for jj=2:(length(z_grid)-1)
%         pi_z(ii,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2-Params.gamma*z_grid(ii),0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2-Params.gamma*z_grid(ii),0,sigma);
%     end
%     pi_z(ii,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2-Params.gamma*z_grid(ii),0,sigma);
% end
% z_grid=gpuArray(z_grid);
% pi_z=gpuArray(pi_z);
% Double check: cumsum(pi_z,2) shows all each row adding to one.

%% General eqm variables: give some initial values
GEPriceParamNames={'r','b','T'};
Params.r=0.06; % interest rate on assets
Params.b=1.2; % Benefits level for retirees
Params.T=1; % lumpsum transfers made out of the accidental bequests
% I originally had b=0.8, T=0.6; switched to these after solving for the GE as I know they are
% closer to the true values, which helps make things run a bit faster.

% Following are coded into the return function to get the values of w and tau given the value of r
% KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% KdivY=(KdivL^(1-Params.alpha))/Params.A;
% Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % wage rate (per effective labour unit)
% Params.tau=0.195*(1-Params.delta*KdivY);

% Note that G is not part of the GEPriceParamNames, this is because it is
% effectively just a residual of the model and plays no part in the actual
% computation. It doesn't effect anyones behaviour, so we don't need it to
% solve the model, and once we solve the model if we want to know what is it
% we can just calculate it based on the governments budget balance equation.


%% Grids
% w is dec in r, so on the assumption that r will not be any lower than zero, we know that w will not be higher than:
maxw=Params.A*(1-Params.alpha)*((((0+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1)))^Params.alpha);
maxa=250;
a_grid=[linspace(-maxw,15,ceil(n_a/2))'; linspace(15,maxa,n_a-ceil(n_a/2))']; % Chose 15 as this is anyway greater than the mean
a_grid(ceil(n_a/2))=0; %There are two 15 values, replace one with 0 
a_grid=sort(a_grid); % Double-check: length(unique(a_grid))==n_a
% Note that the borrowing constraint is instead enforced inside the return
% function. This means some grid points are wasted, but was a bit cleaner.

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};
 
ReturnFn=@(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw) Huggett1996_ReturnFn(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw)
ReturnFnParamNames={'sigma','r','ybarj','theta','b','bvec','T','delta','alpha','A','bc_equalsminusw'}; %It is important that these are in same order as they appear in 'Huggett1996_ReturnFn'

vfoptions.verbose=0;
vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)


%% Agents age distribution
% Almost all OLG models include some kind of population growth, and perhaps
% some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
Params.mewj=ones(1,Params.J);
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj);

simoptions.nsims=4*10^5;
simoptions.iterate=1;
AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

Params.fractionretired=sum(Params.mewj.*Params.bvec); % Note: bvec is really just an indicator of retirement

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Steady State Aggregates (important that ordering of Names and Functions is the same)
SSvaluesParamNames=struct();
SSvaluesParamNames(1).Names={};
SSvaluesFn_1 = @(aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
SSvaluesParamNames(2).Names={'ybarj'};
SSvaluesFn_2 = @(aprime_val,a_val,z_val,ybarj) exp(z_val+ybarj); % Aggregate labour supply (in efficiency units)
SSvaluesParamNames(3).Names={'sj','r','tau'};
SSvaluesFn_3 = @(aprime_val,a_val,z_val,sj,r,tau) (1-sj)*aprime_val*(1+r*(1-tau)); % Total accidental bequests
SSvaluesFn={SSvaluesFn_1,SSvaluesFn_2,SSvaluesFn_3};
% Note that the aggregate labour supply is actually entirely exogenous and so I could just precompute it, but am feeling lazy.

% General Equilibrium Equations
% Recall that GEPriceParamNames={'r','b','T'}; In following lines p is the vector of these and so, e.g., p(2) is T.
GeneralEqmEqnParamNames=struct();
GeneralEqmEqnParamNames(1).Names={'A','alpha','delta'};
GeneralEqmEqn_1 = @(AggVars,p,A,alpha,delta) p(1)-(A*alpha*(AggVars(1)^(alpha-1))*(AggVars(2)^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqnParamNames(2).Names={'theta','fractionretired','A','alpha'};
GeneralEqmEqn_2 = @(AggVars,p,theta,fractionretired,alpha,A) p(2)*fractionretired-theta*(A*(1-alpha)*(AggVars(1)^(alpha))*(AggVars(2)^(-alpha)))*AggVars(2); % Retirement benefits equal Payroll tax revenue: b*fractionretired-theta*w*L
GeneralEqmEqnParamNames(3).Names={};
GeneralEqmEqn_3 = @(AggVars,p) p(3)-AggVars(3); % Lump-sum transfers equal Accidental bequests 
GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_2,GeneralEqmEqn_3};

%% We have finished setting up the model, now we need to solve it for all the variations that Huggett (1996) does.
% Hugget solves more than one economy
% Params.sigma=1.5, 3
sigmavec=[1.5,3];
% Params.alowerbar=0, -w
bc_equalsminuswvec=[0,1];
% Params.sigmasqepsilon=0, 0.045
sigmasqepsilonvec=[0, 0.045];
% Certain lifetimes: beta=0.994, Params.sj=1, beta=1.011, Params.sj=vector
for sigma_c=1:2
    Params.sigma=sigmavec(sigma_c);
    for alowerbar_c=1:2
        Params.bc_equalsminusw=bc_equalsminuswvec(alowerbar_c);
        for sigmasqepsilon_c=1:2
            Params.sigmasqepsilon=sigmasqepsilonvec(sigmasqepsilon_c);
            for uncertainlifetime_c=1:2
                fprintf('Current economy is %i %i %i %i \n', sigma_c, alowerbar_c, sigmasqepsilon_c, uncertainlifetime_c)
                if uncertainlifetime_c==1 % Uncertain lifetime
                    Params.beta=1.011; 
                    Params.sj=1-Params.dj;
                elseif uncertainlifetime_c==2 % Certain lifetime
                    Params.beta=0.994;
                    Params.sj=ones(size(Params.dj));
                end
                OutputResults=Huggett1996_Fn(Params, n_a,n_z,N_j, a_grid, ReturnFn, DiscountFactorParamNames, ReturnFnParamNames, AgeWeightsParamNames, SSvaluesFn, SSvaluesParamNames, GEPriceParamNames, GeneralEqmEqns, GeneralEqmEqnParamNames, simoptions,vfoptions,CheckUniquenessOfGE);
                FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults=OutputResults;
                counter=[sigma_c, alowerbar_c, sigmasqepsilon_c, uncertainlifetime_c];
                save ./SavedOutput/Huggett1996_Counter.mat counter
            end
        end
    end
end

save ./SavedOutput/Huggett1996_FullResults.mat FullResults
% load ./SavedOutput/Huggett1996_FullResults.mat FullResults

%% Draw figures from Huggett (1996)

% Fig 1
figure(1);
plot(19+(1:1:Params.J),exp(Params.ybarj)./sum(exp(Params.ybarj).*Params.mewj))
title({'Earnings Profile (ratio to overall mean)'})
xlabel('Age')
ylabel('Earnings')
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure1.png')

AgeConditionalStats=FullResults(1,2,2,1).OutputResults.AgeConditionalStats;
% Fig 2
figure(2);
plot(19+(1:1:71),AgeConditionalStats.Mean(1:end-8),19+(1:1:71),AgeConditionalStats.QuantileCutoffs(2,1:end-8),19+(1:1:71),AgeConditionalStats.QuantileCutoffs(5,1:end-8),19+(1:1:71),AgeConditionalStats.QuantileCutoffs(10,1:end-8))
legend('Mean','10th','25th','Median')
title({'Wealth Profiles: Uncertain Lifetimes'})
xlabel('Age')
ylabel('Wealth')
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure2.png')

% Figure 3 of Huggett (1996) is based on US data.

% Fig 4
figure(3);
AgeConditionalStats1=FullResults(1,1,2,2).OutputResults.AgeConditionalStats;
AgeConditionalStats2=FullResults(1,2,2,2).OutputResults.AgeConditionalStats;
plot(19+(11:1:56),AgeConditionalStats1.Gini(11:(end-8-15)),19+(11:1:56),AgeConditionalStats2.Gini(11:(end-8-15)))
legend('a=0, sigma=1.5', 'a=-w, sigma=1.5')
title({'Certain Lifetimes: Gini coefficients within age groups'})
xlabel('Age')
ylabel('Wealth Gini')
ylim([0,1.4])
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure4.png')

% Fig 5
figure(4);
AgeConditionalStats1=FullResults(1,1,2,1).OutputResults.AgeConditionalStats;
AgeConditionalStats2=FullResults(1,2,2,1).OutputResults.AgeConditionalStats;
plot(19+(11:1:56),AgeConditionalStats1.Gini(11:(end-8-15)),19+(11:1:61),AgeConditionalStats2.Gini(11:(end-8-10)))
legend('a=0, sigma=1.5', 'a=-w, sigma=1.5')
title({'Uncertain Lifetimes: Gini coefficients within age groups'})
xlabel('Age')
ylabel('Wealth Gini')
ylim([0,1.4])
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure5.png')


%% Tables from Huggett (1996)

Table3=cell(8,9);
sigma_c=1;
alowerbarstr={'0','-w'};
% Row 1
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(1,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 2
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(2,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 3
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(3,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 4
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(4,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 5
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(5,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 6
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(6,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 7
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(7,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 8
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(8,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};

%Table 3
FID = fopen('./SavedOutput/LatexInputs/Huggett1996_Table3.tex', 'w');
fprintf(FID, '\\center{Wealth Distribution (risk aversion coefficient $\\sigma=$%8.1f) \\\\ \n ', sigmavec(sigma_c));
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllllll} \\hline \n');
fprintf(FID, ' &  & &  &  & \\multicolumn{3}{l}{Percentage wealth in} & \\\\ \n');
fprintf(FID, 'Credit & Earnings &       & Transfer &         & \\multicolumn{3}{l}{the top} & Zero or \\\\ \n');
fprintf(FID, 'limit  & shock    &       & wealth   & Wealth  &                              & negative \\\\ \n');
fprintf(FID, '$\\underbar{a}$  & $\\sigma_{\\epsilon}^2$ & $K/Y$ & ration & Gini & 1\\%% & 5\\%% & 20 \\%% & wealth (\\%%) \\\\ \n');
fprintf(FID, '\\multicolumn{2}{l}{US Economy} & 3.0 & 0.78--1.32 & 0.72 & 28 & 49 & 75 & 5.8--15.0 \\\\ \n');
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Certain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{1,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{2,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{3,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{4,:});
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Uncertain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{5,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{6,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{7,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{8,:});
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d. \\\\ \n', n_a, n_z);
fprintf(FID, '}} \\end{minipage} }');
fclose(FID);



Table4=cell(8,9);
sigma_c=2;
alowerbarstr={'0','-w'};
% Row 1
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(1,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 2
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(2,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 3
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(3,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 4
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(4,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 5
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(5,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 6
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(6,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 7
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(7,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 8
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(8,:)={alowerbarstr{alowerbar_c},sigmaepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};

%Table 3
FID = fopen('./SavedOutput/LatexInputs/Huggett1996_Table4.tex', 'w');
fprintf(FID, '\\center{Wealth Distribution (risk aversion coefficient $\\sigma=$%8.1f) \\\\ \n ', sigmavec(sigma_c));
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllllll} \\hline \n');
fprintf(FID, ' &  & &  &  & \\multicolumn{3}{l}{Percentage wealth in} & \\\\ \n');
fprintf(FID, 'Credit & Earnings &       & Transfer &         & \\multicolumn{3}{l}{the top} & Zero or \\\\ \n');
fprintf(FID, 'limit  & shock    &       & wealth   & Wealth  &                              & negative \\\\ \n');
fprintf(FID, '$\\underbar{a}$  & $\\sigma_{\\epsilon}^2$ & $K/Y$ & ration & Gini & 1\\%% & 5\\%% & 20 \\%% & wealth (\\%%) \\\\ \n');
fprintf(FID, '\\multicolumn{2}{l}{US Economy} & 3.0 & 0.78--1.32 & 0.72 & 28 & 49 & 75 & 5.8--15.0 \\\\ \n');
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Certain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{1,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{2,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{3,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{4,:});
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Uncertain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{5,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{6,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{7,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{8,:});
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d. \\\\ \n', n_a, n_z);
fprintf(FID, '}} \\end{minipage} }');
fclose(FID);





























