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



% A few lines I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
try % Server has 16 cores, but is shared with other users, so use max of 8.
    parpool(8)
    gpuDevice(1)
catch % Desktop has less than 8, so will give error, on desktop it is fine to use all available cores.
    parpool
end
PoolDetails=gcp;
NCores=PoolDetails.NumWorkers;

%% State space variables of the model
% Age
% One endogenous variable: assets
% Two stochastic exogenous variables: income shock, medical expense shock
% Three fixed types: 
%  ft1: high-school dropout
%  ft2: high-school
%  ft3: university

%% Declare some parameter values

Params.J=80; % Ages 21 to 100 inclusive.
WorkingAgeVec=21:1:65;
RetiredAgeVec=66:100;
Params.age=[WorkingAgeVec,RetiredAgeVec];

% 
n_a=751;
maxa=5*10^5;
n_z=[15,15]; % income, medical
N_i=3; % Number of fixed types
N_j=Params.J; % Number of periods in finite horizon
Params.q=3; % For tauchen method. They used q=2.5: pg 112 of pdf, "range between 2.5 standard deviations (of the unconditional distribution) above and below..."

Params.gamma=3;
% Gamma plays three roles:      
        % gamma is coefficient of relative risk aversion
        % 1/gamma is the intertemporal elasticity of substitution of consumption
        % gamma+1 is the coefficient of relative prudence (Kimball, 1990)

Params.delta=0.03; % rate at which agents discount future
Params.beta=1/(1+Params.delta);
Params.r=0.03; % interest rate on assets

% Mortality rates (probability of survival)
% Table A.1
Params.dj=[0.00058, 0.00061, 0.00062, 0.00064, 0.00065, 0.00067, 0.00069, 0.00070, 0.00072, 0.00075, 0.00078, 0.00082, 0.00086, 0.00091, 0.00098, 0.00105, 0.00115, 0.00128, 0.00144, 0.00161, 0.00180, 0.00200,...
    0.00221, 0.00242, 0.00266, 0.00292, 0.00320, 0.00349, 0.00380, 0.00413, 0.00450, 0.00490, 0.00533,...
    0.00533, 0.00581, 0.00632, 0.00689, 0.00749, 0.00811, 0.00878, 0.00952, 0.01033, 0.01121, 0.01223,...
    0.01332, 0.01455, 0.01590, 0.01730, 0.01874, 0.02028, 0.02203, 0.02404, 0.02623, 0.02863, 0.03128,...
    0.03432, 0.03778, 0.04166, 0.04597, 0.05078, 0.05615, 0.06214, 0.06885, 0.07631, 0.08455, 0.09352,...
    0.10323, 0.11367, 0.12484, 0.13677, 0.14938, 0.16289, 0.17721, 0.19234, 0.20828, 0.22418, 0.23980, 0.25495, 0.26937, 0.28284]; % conditional probability of death
Params.sj=1-Params.dj; % Conditional survival probabilities. Act as a discount rate.

% Wage (Earnings) incomes, W
% Table A.2
Params.DeterministicWj_working.ft1=[10993-196*WorkingAgeVec+2167*((WorkingAgeVec.^2)/100)-2655*((WorkingAgeVec.^3)/10000)];
Params.DeterministicWj_working.ft2=[-11833+1421*WorkingAgeVec-137*((WorkingAgeVec.^2)/100)-2186*((WorkingAgeVec.^3)/10000)]; 
Params.DeterministicWj_working.ft3=[72270-5579*WorkingAgeVec+19200*((WorkingAgeVec.^2)/100)-18076*((WorkingAgeVec.^3)/10000)];
% Table A.3
Params.DeterministicWj_retired.ft1=[17861-133*RetiredAgeVec]; 
Params.DeterministicWj_retired.ft2=[29733-245*RetiredAgeVec]; 
Params.DeterministicWj_retired.ft3=[48123-429*RetiredAgeVec];
% Is not mentioned in paper, but based on Figures 3a-c is clear that
% earnings from ages 90+ are simply frozen at their age 89 level.
Params.DeterministicWj_retired.ft1(25:35)=Params.DeterministicWj_retired.ft1(24)*ones(size(Params.DeterministicWj_retired.ft1(25:35)));
Params.DeterministicWj_retired.ft2(25:35)=Params.DeterministicWj_retired.ft2(24)*ones(size(Params.DeterministicWj_retired.ft2(25:35)));
Params.DeterministicWj_retired.ft3(25:35)=Params.DeterministicWj_retired.ft3(24)*ones(size(Params.DeterministicWj_retired.ft3(25:35)));
Params.DeterministicWj.ft1=[Params.DeterministicWj_working.ft1, Params.DeterministicWj_retired.ft1];
Params.DeterministicWj.ft2=[Params.DeterministicWj_working.ft2, Params.DeterministicWj_retired.ft2];
Params.DeterministicWj.ft3=[Params.DeterministicWj_working.ft3, Params.DeterministicWj_retired.ft3];
% Compare to Fig A.1(a-c): they won't be exactly the same as drop the year fixed effects, but eyeballing suggests they are fine.
plot(21:1:100, Params.DeterministicWj.ft1, 21:1:100, Params.DeterministicWj.ft2, 21:1:100, Params.DeterministicWj.ft3)
% Stochastic Wj (Table A.4): u_it, an AR(1) in regressions in paper
Params.w_rho=[0.955; 0.946; 0.955];
Params.w_sigmasqepsilon=[0.033; 0.025; 0.016];
Params.w_sigmasqu=Params.w_sigmasqepsilon./(1-Params.w_rho.^2);
Params.w_sigmasqupsilon=[0.040; 0.021; 0.014]; % Estimated from PSID but not used in model.
% Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% "This specification ensures that when we compare the certainty case with
% the earnings uncertainty case, we hold the age-conditional means of earnings constant."

% Medical Expenses, M
% Table A.5
Params.DeterministicMj_working=[5.373+0.073*WorkingAgeVec-0.753*((WorkingAgeVec.^2)/1000); 4.749+0.106*WorkingAgeVec-1.084*((WorkingAgeVec.^2)/1000); 4.816+0.109*WorkingAgeVec-1.090*((WorkingAgeVec.^2)/1000)];
Params.DeterministicMj_retired=[14.441-0.200*RetiredAgeVec+1.288*((RetiredAgeVec.^2)/1000); 11.371-0.101*RetiredAgeVec+0.540*((RetiredAgeVec.^2)/1000); 9.553-0.054*RetiredAgeVec-0.297*((RetiredAgeVec.^2)/1000)];
Params.DeterministicMj=[Params.DeterministicMj_working, Params.DeterministicMj_retired]; % is in logs
Params.m_rho=0.901*ones(3,1);
Params.m_sigmasqepsilon=[0.175; 0.156; 0.153];
Params.m_sigmasqmew=Params.m_sigmasqepsilon./(1-Params.m_rho.^2);
Params.m_sigmasqomega=[0.220; 0.220; 0.220];
% Consumption Floor
Params.Cbar=7000; % (middle of pg. 111)

%% Create some outputs to replicate those of HubbardSkinnerZeldes1994

% Note: The annual growth rate of the population is assumed to be 1%.
Params.n=0.01;

simoptions.numbersims=10^5; % Number of simulations on which panel data (and life-cycle profile) results will be based.

PTypeDist=[0.22,0.56,0.22]'; % Hubbard, Skinner & Zeldes (1994) do not appear to report these 
                             % weights/probabilities anywhere in either of the two papers. 
                             % I have 'ballparked' them based on alternative sources for 1984 as
                             % fractions of US population (pg 1 of http://www.russellsage.org/sites/all/files/chartbook/Educational%20Attainment%20and%20Achievement.pdf )


%% Create the replication results
% Note that we here solve for many parameter/shock cases for which no
% results are actually reported in paper.
Table1=struct();
Table2=nan(12,6);

temp_n_z=n_z; % Store this so that can restore it later after looking at 'certainty cases'

% Get results for uncertainty model
deltavec=[0.03,0.1,0.015];
gammavec=[1,3,5];
Cbarvec=[1,7000];
for delta_c=1:length(deltavec)
    for gamma_c=1:length(gammavec)
        for Cbar_c=1:2
            fprintf('Currently solving for delta_c=%d, gamma_c=%d,Cbar_c=%d',delta_c,gamma_c,Cbar_c)
            Params.Cbar=Cbarvec(Cbar_c);
            Params.delta=deltavec(delta_c);
            Params.gamma=gammavec(gamma_c);
            descriptivestr=['Cbar',num2str(Params.Cbar),'gamma',num2str(Params.gamma),'delta',num2str(Params.delta)];
            descriptivestr(descriptivestr=='.') = []; % Get rid of decimal points
            descriptivestr={descriptivestr};
            [Table1row.(descriptivestr{:}), Table2temp, Table3temp, LifeCycProfiles.(descriptivestr{:})]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
            if Params.Cbar==7000 && Params.delta==0.03 && Params.gamma==3
                Table2(:,1:3)=Table2temp;
            elseif Params.Cbar==1 && Params.delta==0.1 && Params.gamma==3
                Table2(:,4:6)=Table2temp;
            end
            if Params.Cbar==7000 && Params.delta==0.03 && Params.gamma==3
                Table3=Table3temp;
            end
        end
    end
end

% Set Parameters back to defaults
Params.delta=0.03;
Params.gamma=3;
Params.Cbar=7000;

% Get results for 'only lifetime uncertain' case
n_z=[1,1];  Params.w_sigmasqu=[0;0;0];  Params.m_sigmasqmew=[0;0;0];
Params.Cbar=7000;
% [Table1row.allcertain, ~, ~, LifeCycProfiles.(descriptivestr{:})]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
for delta_c=1:length(deltavec)
    for gamma_c=1:length(gammavec)
        for Cbar_c=1:2
            Params.Cbar=Cbarvec(Cbar_c);
            Params.delta=deltavec(delta_c);
            Params.gamma=gammavec(gamma_c);
            descriptivestr=['onlylifetimeuncertain_','Cbar',num2str(Params.Cbar),'gamma',num2str(Params.gamma),'delta',num2str(Params.delta)];
            descriptivestr(descriptivestr=='.') = []; % Get rid of decimal points
            descriptivestr={descriptivestr};
            [Table1row.(descriptivestr{:}), ~, ~, LifeCycProfiles.(descriptivestr{:})]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
        end
    end
end

% Get results for 'certain lifetimes' case
Params.sj=ones(size(Params.sj));
% Paper doesn't mention it but based on Figures it is clear that under 'certain lifetimes' death occours at age 80.
Params.sj((end-19):end)=0;
n_z=temp_n_z; Params.w_sigmasqu=Params.w_sigmasqepsilon./(1-Params.w_rho.^2);  Params.m_sigmasqmew=Params.m_sigmasqepsilon./(1-Params.m_rho.^2);
Params.Cbar=1;
% [Table1row.certainlifetimes, ~, ~, LifeCycProfiles.certainlifetimes]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
for delta_c=1:length(deltavec)
    for gamma_c=1:length(gammavec)
        Params.delta=deltavec(delta_c);
        Params.gamma=gammavec(gamma_c);
        descriptivestr=['certainlifetimes_','Cbar',num2str(Params.Cbar),'gamma',num2str(Params.gamma),'delta',num2str(Params.delta)];
        descriptivestr(descriptivestr=='.') = []; % Get rid of decimal points
        descriptivestr={descriptivestr};
        [Table1row.(descriptivestr{:}), ~, ~, LifeCycProfiles.(descriptivestr{:})]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
    end
end

% Get results for 'all certain' case
Params.sj=ones(size(Params.sj)); 
% Paper doesn't mention it but based on Figures it is clear that under 'certain lifetimes' death occours at age 80.
Params.sj((end-19):end)=0;
n_z=[1,1]; Params.w_sigmasqu=[0;0;0]; Params.m_sigmasqmew=[0;0;0];
Params.Cbar=1;
% [Table1row.allcertain, ~, ~, LifeCycProfiles.(descriptivestr{:})]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
for delta_c=1:length(deltavec)
    for gamma_c=1:length(gammavec)
        Params.delta=deltavec(delta_c);
        Params.gamma=gammavec(gamma_c);
        descriptivestr=['allcertain_','Cbar',num2str(Params.Cbar),'gamma',num2str(Params.gamma),'delta',num2str(Params.delta)];
        descriptivestr(descriptivestr=='.') = []; % Get rid of decimal points
        descriptivestr={descriptivestr};
        [Table1row.(descriptivestr{:}), ~, ~, LifeCycProfiles.(descriptivestr{:})]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
    end
end

n_z=temp_n_z; % Restore this now that we are done looking at all the certainty cases.

save ./SavedOutput/HSZ1994_Tables.mat Table1row Table2 Table3 LifeCycProfiles Params

%% Print out results for the three Tables
load ./SavedOutput/HSZ1994_Tables.mat Table1row Table2 Table3 LifeCycProfiles Params

%Table 1
FID = fopen('./SavedOutput/HubbardSkinnerZeldes_Table1.tex', 'w');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccccc} \n \\hline \\hline \n');
fprintf(FID, ' & &    & \\multicolumn{4}{c}{----Asset-Income Ratio----} & $\\quad$ & \\multicolumn{4}{c}{----Savings Rate----}\\\\ \n');
fprintf(FID, ' & \\multicolumn{2}{c}{Parameter}   & No High & High  &         & Aggre- & $\\quad$ & No High & High &          & Aggre- \\\\ \n');
fprintf(FID, ' & \\multicolumn{2}{c}{Assumptions} & School & School & College & gate & $\\quad$ & School & School & College & gate \\\\ \n');
fprintf(FID, 'Certain Lifespan,   & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.certainlifetimes_Cbar1gamma3delta003);
fprintf(FID, 'earnings, and       & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.certainlifetimes_Cbar1gamma1delta003);
fprintf(FID, 'out of pocket       & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.certainlifetimes_Cbar1gamma5delta003);
fprintf(FID, 'medical             & $\\delta=%8.3f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.certainlifetimes_Cbar1gamma3delta0015);
fprintf(FID, 'expenses            & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.certainlifetimes_Cbar1gamma3delta01);
fprintf(FID, '                    &  & & & & & & & & & & \\\\ \n');
fprintf(FID, 'Uncertain Lifespan, & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.Cbar1gamma3delta003);
fprintf(FID, 'earnings, and       &  & & & & & & & & & & \\\\ \n');
fprintf(FID, 'medical             &  & & & & & & & & & & \\\\ \n');
fprintf(FID, 'expenses            &  & & & & & & & & & & \\\\ \n');
fprintf(FID, '$\\bar{C}=\\$1$     &  & & & & & & & & & & \\\\ \n');
fprintf(FID, '                    &  & & & & & & & & & & \\\\ \n');
fprintf(FID, 'Uncertain Lifespan, & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.Cbar7000gamma3delta003);
fprintf(FID, 'earnings, and       & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.Cbar7000gamma1delta003);
fprintf(FID, 'medical             & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.Cbar7000gamma5delta003);
fprintf(FID, 'expenses            & $\\delta=%8.3f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.Cbar7000gamma3delta0015);
fprintf(FID, '$\\bar{C}=\\$7000$  & $\\delta=%8.2f$, & $\\gamma=%d$ & %8.2f & %8.2f & %8.2f & %8.2f & & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table1row.Cbar7000gamma3delta01);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d and %d. \\\\ \n', n_a, n_z(1), n_z(2));
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%Table 2
FID = fopen('./SavedOutput/HubbardSkinnerZeldes_Table2.tex', 'w');
fprintf(FID, '\\center{Percentage of Households with Consumption Approximately Equal to Income \\\\ \n (Absolute Average Savings Rate $<$ 0.5 Percent of Income} \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{3}{c}{PSID} & \\multicolumn{3}{c}{Simulated}  & \\multicolumn{3}{c}{Simulated} \\\\ \n');
fprintf(FID, ' & \\multicolumn{3}{c}{}     & \\multicolumn{3}{c}{$\\delta=0.03$, Floor=\\$7000}  & \\multicolumn{3}{c}{$\\delta=0.10$, Floor=\\$1} \\\\ \n');
fprintf(FID, 'Age   & NHS & HS & Col. & NHS & HS & Col. & NHS & HS & Col. \\\\ \n');
fprintf(FID, '$<$29   & 0.362 & 0.059 & 0.060 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(1,:));
fprintf(FID, '30-39 & 0.157 & 0.064 & 0.040 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(2,:));
fprintf(FID, '40-49 & 0.067 & 0.017 & 0.025 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(3,:));
fprintf(FID, '50-59 & 0.103 & 0.032 & 0.011 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(4,:));
fprintf(FID, '60-69 & 0.095 & 0.020 & 0.038 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(5,:));
fprintf(FID, 'Total & 0.116 & 0.038 & 0.031 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(6,:));
fprintf(FID, ' & \\multicolumn{9}{c}{For Households with Initial Assets $<$ 0.5 x Average Income}  \\\\ \n');
fprintf(FID, '$<$29   & 0.389 & 0.080 & 0.069 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(7,:));
fprintf(FID, '30-39 & 0.205 & 0.101 & 0.052 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(8,:));
fprintf(FID, '40-49 & 0.133 & 0.045 & 0.031 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(9,:));
fprintf(FID, '50-59 & 0.272 & 0.114 & 0.135 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(10,:));
fprintf(FID, '60-69 & 0.302 & 0.083 & 0.381 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(11,:));
fprintf(FID, 'Total & 0.252 & 0.087 & 0.060 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table2(12,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'NHS=No high-school, HS=High-school, Col.=College. Numbers for PSID are those of original study, not part of replication. \\\\ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d and %d. \\\\ \n', n_a, n_z(1), n_z(2));
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 3
FID = fopen('./SavedOutput/HubbardSkinnerZeldes_Table3.tex', 'w');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}ccccc} \n \\hline \\hline \n');
fprintf(FID, ' Dependent Varialbe & $\\Delta C*$ & $\\Delta ln(C)*$ & $\\Delta ln(C)**$ & $\\Delta ln(C)**$\\\\ \n');
fprintf(FID, ' & \\multicolumn{4}{c}{$\\delta=0.03$, Floor=\\$7000} \\\\ \n');
fprintf(FID, '$\\Delta Y$     & %8.3f   & & & \\\\ \n', Table3(1,1));
fprintf(FID, '                & (%8.2f) & & & \\\\ \n', Table3(2,1));
fprintf(FID, '$\\Delta ln(Y)$ & & %8.3f & %8.3f & %8.3f \\\\ \n', Table3(3,2:4));
fprintf(FID, '                & & (%8.2f) & (%8.2f) & (%8.2f) \\\\ \n', Table3(4,2:4));
fprintf(FID, '$Age$           & & & & %8.3f \\\\ \n', Table3(5,4));
fprintf(FID, '                & & & & (%8.2f) \\\\ \n', Table3(6,4));
fprintf(FID, '$Age^2$         & & & & %8.3f  \\\\ \n', Table3(7,4));
fprintf(FID, '                & & & & (%8.2f)  \\\\ \n', Table3(8,4));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Source: Simulated data from model under the benchmark case ($\\delta = 0.03$, $\\gamma = 3$) with a consumption floor of \\$7000. Absolute values of t-statistics are in parentheses. The Campbell and Mankiw (1989) coefficient (in levels, corresponding to the first column) is 0.469, and the Lusardi (1993) coefficient (in logs, corresponding to columns 2 through 4) is 0.409. \\\\ \n');
fprintf(FID, '*: Instruments are two and three year lags of consumption and income, as well as age and age-squared. \\\\ \n');
fprintf(FID, '**: Instruments are one, two and three year lags of consumption and income, as well as age and age-squared. \\\\ \n');
fprintf(FID, 'Note: Original regressions by Cambell-Mankiw were on aggregate data. Here are on microdata. \\\\ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d and %d. \\\\ \n', n_a, n_z(1), n_z(2));
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Generate the Figures

% Figure 1
figure(1)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft1(1,1:(end-19),1)/1000)
hold on
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft2(1,1:(end-19),1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft3(1,1:(end-19),1)/1000)
hold off
title({'Average Assets by Age';'All Certain, $1 Floor'})
legend('No High School Degree', 'High School Degree', 'College')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure1.png')

% Figure 2a
figure(2)
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft1(1,:,1)/1000)
hold on
plot(Params.age, LifeCycProfiles.Cbar7000gamma3delta003.ft1(1,:,1)/1000)
plot(Params.age, LifeCycProfiles.onlylifetimeuncertain_Cbar1gamma3delta003.ft1(1,:,1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft1(1,1:(end-19),1)/1000)
hold off
title({'Average Assets by Age';'No High School Degree'})
legend('All uncertain ($1 floor)', 'All uncertain ($7000 floor)','Only lifetime uncertain', 'All certain')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure2a.png')


% Figure 2b
figure(3)
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft2(1,:,1)/1000)
hold on
plot(Params.age, LifeCycProfiles.Cbar7000gamma3delta003.ft2(1,:,1)/1000)
plot(Params.age, LifeCycProfiles.onlylifetimeuncertain_Cbar1gamma3delta003.ft2(1,:,1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft2(1,1:(end-19),1)/1000)
hold off
title({'Average Assets by Age';'High School Degree'})
legend('All uncertain ($1 floor)', 'All uncertain ($7000 floor)','Only lifetime uncertain', 'All certain')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure2b.png')


% Figure 2c
figure(4)
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft3(1,:,1)/1000)
hold on
plot(Params.age, LifeCycProfiles.Cbar7000gamma3delta003.ft3(1,:,1)/1000)
plot(Params.age, LifeCycProfiles.onlylifetimeuncertain_Cbar1gamma3delta003.ft3(1,:,1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft3(1,1:(end-19),1)/1000)
hold off
title({'Average Assets by Age';'College Degree'})
legend('All uncertain ($1 floor)', 'All uncertain ($7000 floor)','Only lifetime uncertain', 'All certain')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure2c.png')



% Figure 3a
figure(5)
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft1(2,:,1)/1000) % Earnings
hold on
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft1(2,1:(end-19),1)/1000) % Earnings
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft1(5,:,1)/1000) % Consumption
plot(Params.age, LifeCycProfiles.Cbar7000gamma3delta003.ft1(5,:,1)/1000)
plot(Params.age, LifeCycProfiles.onlylifetimeuncertain_Cbar1gamma3delta003.ft1(5,:,1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft1(5,1:(end-19),1)/1000)
hold off
title({'Average Consumption and Earnings by Age';'No High School Degree'})
legend('Earnings: All Uncertain','Earnings: All certain','Cons: All uncertain ($1 floor)', 'Cons: All uncertain ($7000 floor)','Cons: Only lifetime uncertain', 'Cons: All certain')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure3a.png')


% Figure 3b
figure(6)
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft2(2,:,1)/1000) % Earnings
hold on
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft2(2,1:(end-19),1)/1000) % Earnings
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft2(5,:,1)/1000) % Consumption
plot(Params.age, LifeCycProfiles.Cbar7000gamma3delta003.ft2(5,:,1)/1000)
plot(Params.age, LifeCycProfiles.onlylifetimeuncertain_Cbar1gamma3delta003.ft2(5,:,1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft2(5,1:(end-19),1)/1000)
hold off
title({'Average Consumption and Earnings by Age';'High School Degree'})
legend('Earnings: All Uncertain','Earnings: All certain','Cons: All uncertain ($1 floor)', 'Cons: All uncertain ($7000 floor)','Cons: Only lifetime uncertain', 'Cons: All certain')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure3b.png')


% Figure 3c
figure(7)
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft3(2,:,1)/1000) % Earnings
hold on
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft3(2,1:(end-19),1)/1000) % Earnings
plot(Params.age, LifeCycProfiles.Cbar1gamma3delta003.ft3(5,:,1)/1000) % Consumption
plot(Params.age, LifeCycProfiles.Cbar7000gamma3delta003.ft3(5,:,1)/1000)
plot(Params.age, LifeCycProfiles.onlylifetimeuncertain_Cbar1gamma3delta003.ft3(5,:,1)/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.allcertain_Cbar1gamma3delta003.ft3(5,1:(end-19),1)/1000)
hold off
title({'Average Consumption and Earnings by Age';'College Degree'})
legend('Earnings: All Uncertain','Earnings: All certain','Cons: All uncertain ($1 floor)', 'Cons: All uncertain ($7000 floor)','Cons: Only lifetime uncertain', 'Cons: All certain')
xlabel('Age')
ylabel('Thousands')
saveas(gcf,'./SavedOutput/Graphs/HubbardSkinnerZeldes1994_Figure3c.png')




