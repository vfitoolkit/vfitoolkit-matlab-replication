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
n_a=250;
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
Params.m_sigmasqomega=[0.220; 0.220; 0.220];
% Consumption Floor
Params.Cbar=7000; % (middle of pg. 111)

%% Create some outputs to replicate those of HubbardSkinnerZeldes1994

% Note: The annual growth rate of the population is assumed to be 1%.
Params.n=0.01;

simoptions.numbersims=10^3; % Number of simulations on which panel data (and life-cycle profile) results will be based.

PTypeDist=[0.22,0.56,0.22]'; % Hubbard, Skinner & Zeldes (1994) do not appear to report these 
                             % weights/probabilities anywhere in either of the two papers. 
                             % I have 'ballparked' them based on alternative sources for 1984 as
                             % fractions of US population (pg 1 of http://www.russellsage.org/sites/all/files/chartbook/Educational%20Attainment%20and%20Achievement.pdf )


%% Create the replication results
Table1=struct();
Table2=nan(12,6);

% Get results for uncertainty model
deltavec=[0.03,0.1,0.15];
gammavec=[1,3,5];
Cbarvec=[1,7000];
for delta_c=1:length(deltavec)
    for gamma_c=1:length(gammavec)
        for Cbar_c=1:2
            fprintf('Currently solving for delta_c=%d, gamma_c=%d,Cbar_c=%d',delta_c,gamma_c,Cbar_c)
            Params.Cbar=Cbarvec(Cbar_c);
            Params.delta=deltavec(delta_c);
            Params.gamma=gammavec(gamma_c);
            descriptivestr={['Cbar',num2str(Params.Cbar),'gammma',num2str(Params.gamma),'delta',num2str(100*Params.delta)]};
            [Table1row.(descriptivestr{:}), Table2temp, Table3temp, LifeCycProfiles]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);
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

% Get results for 'certain lifetimes' case
Params.sj=ones(size(Params.sj));
Params.Cbar=1;
[Table1row.certainlifetimes, ~, ~, LifeCycProfiles.certainlifetimes]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);

% Get results for 'all certainty' case
Params.sj=ones(size(Params.sj)); 
n_z=[1,1];
Params.Cbar=1;
[Table1row.allcertain, ~, ~, LifeCycProfiles.allcertain]=HubbardSkinnerZeldes1994_function(Params,n_a,n_z,simoptions, PTypeDist);

