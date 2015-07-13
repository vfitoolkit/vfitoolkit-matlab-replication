%%% Example using the model of Hansen (1985)
% Note: Most of the run time is actually spent in creating the
% discretization of the shock process [AR(1) shock with log-normal
% innovations means standard methods cannot be used].
n_d=51;  % decide whether to work
n_a=501; % assets
n_z=31;  % tech shock
AlternativeProductivityShocks=0


%% Setup for server
addpath(genpath('./MatlabToolkits/'))
% try % Server has 16 cores, but is shared, so use max of 8.
%     PoolDetails=parpool(8)
% catch % Desktop has less than 8, so will give error, on desktop it is fine to use all available cores.
%     PoolDetails=parpool
% end
% %PoolDetails=gcp;
% NCores=PoolDetails.NumWorkers;

%% Some Toolkit options
Parallel=2; % Use GPU

vfoptions.parallel=Parallel;
simoptions.parallel=Parallel;

%% Setup

%Discounting rate
beta = 0.99;

%Parameter values
alpha = 0.36; % alpha (Hansen refers to this as theta)
A=2;
h0=0.53;
B=-A*log(1-h0)/h0;
rho = 0.95; % rho (Hansen refers to this as gamma)
delta = 0.025; % delta
sigma_epsilon=0.00712;
sigmasq_epsilon=sigma_epsilon^2;
%Step 1: Compute the steady state
K_ss=(alpha/(1/beta-1+delta))^(1/(1-alpha))*2/3; % Hansen footnote 15.
% K_ss=((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));


%Create grids (it is very important that each grid is defined as a column vector)

if AlternativeProductivityShocks==0
    % Creating a grid (and transition matrix) for z is slightly awkward due to epsilon being
    % distributed as lognormal (so cannot just use the Tauchen quadrature like
    % we usually might)
    % Instead, simulate z, then do a quadrature with states chosen by
    % matlab's histcounts. Transition matrix is just done by counting the
    % transitions and then normalizing.
    burnin=10^4;
    T=10^8;
    z_sim=zeros(T,1);
    % Calculate the mean and std dev of the normal shock the exponential of
    % which will be a lognormal with mean 1-rho and std dev sigma_epsilon
    sigma_ln=sqrt(log(1+(sigma_epsilon^2)/((1-rho)^2)));
    mu_ln=log(1-rho)-(sigma_ln^2)/2;
    % check=exp(mu_ln+sigma_ln*randn(10^4,1)); mean(check); std(check); % Works
    temp=1-rho; %The unconditional mean of z
    for tt=1:burnin
        temp=rho*temp+exp(mu_ln+sigma_ln*randn(1,1));
    end
    z_sim(1)=temp;
    for tt=2:T
        z_sim(tt)=rho*z_sim(tt-1)+exp(mu_ln+sigma_ln*randn(1,1));
    end
    [N,edges,bin] = histcounts(z_sim,n_z);
    z_grid=zeros(n_z,1);
    for ii=1:n_z
        z_grid(ii)=mean(z_sim(bin==ii));
    end
    pi_z=zeros(n_z,n_z);
    for tt=2:T
        pi_z(bin(tt-1),bin(tt))=pi_z(bin(tt-1),bin(tt))+1;
    end
    pi_z=pi_z./(sum(pi_z,2)*ones(1,n_z));
elseif AlternativeProductivityShocks==1
    % Use productivity shocks modelled as log(z) is AR(1) with normal dist innovations
    q=3;
    tauchenoptions.parallel=Parallel;
    mcmomentsoptions.parallel=Parallel;
    [z_grid, pi_z]=TauchenMethod(0,sigmasq_epsilon,rho,n_z,q,tauchenoptions); % the AR(1) on log(z)
    z_grid=exp(z_grid); % so z is just exp(log(z))
    %Normalize z by dividing it by Expectation_z (so that the resulting process has expectaion equal to 1.
    [Expectation_z,~,~,~]=MarkovChainMoments(z_grid,pi_z,mcmomentsoptions);
    z_grid=z_grid./Expectation_z;
    [Expectation_z,~,~,~]=MarkovChainMoments(z_grid,pi_z,mcmomentsoptions);
end

a_grid=sort([linspace(0+0.0001,K_ss-0.0001,n_a-floor(n_a/2)),linspace(0,2*K_ss,floor(n_a/2))])';
d_grid=linspace(0,1,n_d)';

z_grid=gpuArray(z_grid);
a_grid=gpuArray(a_grid);
d_grid=gpuArray(d_grid);
pi_z=gpuArray(pi_z);

% n_z=length(z_grid);
% n_a=length(a_grid);
% n_d=length(d_grid);

ReturnFn=@(d_val, aprime_val, a_val, z_val, alpha, delta, A, B, Economy) Hansen1985_ReturnFn(d_val, aprime_val, a_val, z_val, alpha, delta, A, B, Economy);
ReturnFnParams=[alpha, delta, A, B, 0]; % The zero is later replaced with Economy.

%% Solve Model and Generate Table 1 of Hansen (1985)

% Hansen footnote (b) of Table 1 says statistics are based on 100 simulations of 115 periods each.
NSims=100;  %Number of simulations
simoptions.simperiods=115;
StdBusCycleStats=zeros(2,7,2,NSims,'gpuArray');

for Economy=1:2 % Divisible and Indivisible labour respectively
    ReturnFnParams(5)=Economy;
    %% Solve
    disp('Solve value fn problem')
    V0=ones([n_a,n_z],'gpuArray'); %(a,z)
    [V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, beta, ReturnFn,vfoptions,ReturnFnParams);

    disp('Sim time series')

    %No need for asyptotic distribution.
    %SteadyStateDist=SteadyState_Case1_Simulation(Policy,n_d,n_a,n_z,pi_z, simoptions);
    %SteadyStateDist=SteadyState_Case1(SteadyStateDist,Policy,n_d,n_a,n_z,pi_z,simoptions);
    %plot(1:1:n_z, N/sum(N), 1:1:n_z, sum(SteadyStateDist,1)) % Can see that the z discretization is working as two lines are identical.
    
    if AlternativeProductivityShocks==0
        if Economy==1
            save ./SavedOutput/Hansen1985_Economy1.mat V Policy
        elseif Economy==2
            save ./SavedOutput/Hansen1985_Economy2.mat V Policy
        end
    end
    %% Generate Table 1 of Hansen (1985)
    
    for ii=1:NSims
        TimeSeriesIndexes=SimTimeSeriesIndexes_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
        
        %Define the functions which we wish to create time series for (from the TimeSeriesIndexes)
        TimeSeriesFn_1 = @(d_val,aprime_val,a_val,z_val) z_val*(a_val^alpha)*(d_val^(1-alpha)); %Output (from eqn 1)
        TimeSeriesFn_2 = @(d_val,aprime_val,a_val,z_val) z_val*(a_val^alpha)*(d_val^(1-alpha)) -(aprime_val-(1-delta)*a_val); %Consumption (=output-investment, from eqn 2; this formula is valid for both divisible and indivisible labour)
        TimeSeriesFn_3 = @(d_val,aprime_val,a_val,z_val) aprime_val-(1-delta)*a_val; %Investment (from eqn 3)
        TimeSeriesFn_4 = @(d_val,aprime_val,a_val,z_val) a_val; %Capital Stock
        TimeSeriesFn_5 = @(d_val,aprime_val,a_val,z_val) d_val; %Hours
        TimeSeriesFn_6 = @(d_val,aprime_val,a_val,z_val) (z_val*(a_val^alpha)*(d_val^(1-alpha)))/d_val; %Productivity (measured in data as output divided by hours)
        TimeSeriesFn_7 = @(d_val,aprime_val,a_val,z_val) z_val; %Tech Shock (Hansen 1985 does not report this, just for interest)
        
        
        TimeSeriesFn={TimeSeriesFn_1, TimeSeriesFn_2, TimeSeriesFn_3, TimeSeriesFn_4, TimeSeriesFn_5, TimeSeriesFn_6, TimeSeriesFn_7};
        
        TimeSeries=TimeSeries_Case1(TimeSeriesIndexes,Policy, TimeSeriesFn, n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);
        
        [OutputTrend,OutputCyc]=hpfilter(gather(log(TimeSeries(1,:))),1600); % hpfilter() does not yet exist for gpu
        for jj=1:7
            [TimeSeriesTrend,TimeSeriesCyc]=hpfilter(gather(log(TimeSeries(jj,:))),1600); % hpfilter() does not yet exist for gpu
            temp=cov(100*OutputCyc,100*TimeSeriesCyc); % To match treatment of data the model output must be logged and then hpfiltered, then multiply by 100 to express as percent
            StdBusCycleStats(1,jj,Economy,ii)=sqrt(temp(2,2)); % Std dev (column a)
            StdBusCycleStats(2,jj,Economy,ii)=temp(1,2)/(sqrt(temp(1,1))*sqrt(temp(2,2))); % Correlation with output (column b)
        end
    end
    StdBusCycleStats_means=mean(StdBusCycleStats,4); % multiply by 100 to make them percents
    StdBusCycleStats_stddev=std(StdBusCycleStats,0,4);
    
end


%% Print out a file containing Table 1
HansenDataResults=[1.76, 1.00; 1.29, 0.85; 8.60, 0.92; 0.63, 0.04; 1.66, 0.76; 1.18, 0.42];
VarNamesStrVec={'Output', 'Consumption','Investment','Capital Stock','Hours','Productivity'};
if AlternativeProductivityShocks==0
    FID = fopen('./SavedOutput/LatexInputs/Table_Hansen1985.tex', 'w');
elseif AlternativeProductivityShocks==1
    FID = fopen('./SavedOutput/LatexInputs/Table_Hansen1985_AltProdShocks.tex', 'w');
end
fprintf(FID, 'Standard deviations in percent (a) and correlations with output (b) for US and artificial economies. \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \\hline \\hline \n');
fprintf(FID, '           & \\multicolumn{2}{c}{Quarterly U.S. Time Series\\textsuperscript{a}} & \\multicolumn{2}{c}{Economy with} & \\multicolumn{2}{c}{Economy with} \\\\ \n');
fprintf(FID, '           & \\multicolumn{2}{c}{1955:Q3-1984:Q1} & \\multicolumn{2}{c}{divisible labor\\textsuperscript{b}} & \\multicolumn{2}{c}{indivisible labor\\textsuperscript{b}} \\\\ \n');
fprintf(FID, 'Series & (a) & (b) & (a) & (b) & (a) & (b) \\\\ \\hline \n');
for jj=1:6
    VarNameStr=VarNamesStrVec{jj};
    fprintf(FID, [VarNameStr,' & %1.2f & %1.2f & %1.2f (%1.2f) & %1.2f (%1.2f) & %1.2f (%1.2f) & %1.2f (%1.2f) \\\\ \n'], HansenDataResults(jj,1), HansenDataResults(jj,2), StdBusCycleStats_means(1,jj,1), StdBusCycleStats_stddev(1,jj,1), StdBusCycleStats_means(2,jj,1), StdBusCycleStats_stddev(2,jj,1), StdBusCycleStats_means(1,jj,2), StdBusCycleStats_stddev(1,jj,2), StdBusCycleStats_means(2,jj,2), StdBusCycleStats_stddev(2,jj,2) );
end
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
%fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
%fprintf(FID, ['\\textsuperscript{a} The US time series used are real GNP, total consumption expenditures, and gross private domestic investment (all in 1972 dollars). The capital stock series includes non-residential equipment and structures. The hours series includes total hours for persons at work in non-agricultural industries as derived from the \\textit{Current Population Survey}. Productivity is output divided by hours. All series are seasonally adjusted, logged and detrended. \\\\ \n']);
%fprintf(FID, ['\\textsuperscript{b} The standard deviations and correlations with output are sample means of statistics computed for each of 100 simulations. Each simulation consists of 115 periods, which is the same number of periods as the US sample. The numbers in parentheses are sample standard deviations of these statistics. Before computing any statistics each simulated time series was logged and detrended using the same procedure used for the US time series. \\\\ \n']);
%fprintf(FID, '}} \\end{minipage}');
fclose(FID);









