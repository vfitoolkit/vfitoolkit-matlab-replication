% This code replicates the results of Hubbard, Skinner & Zeldes (1994WP) - The importance of precautionary motives in explaining individual and aggregate saving

% Age
% One endogenous variable: assets
% Two stochastic exogenous variables: income shock, medical expense shock

J=80; % Ages 21 to 100 inclusive.
WorkingAgeVec=21:1:65;
RetiredAgeVec=66:100;

% 
n_a=250;
maxa=5*10^5;
n_z=[15,15]; % income, medical
N_i=3; % Number of fixed types
N_j=J; % Number of periods in finite horizon
Params.q=3; % For tauchen method

Params.gamma=3;
% Gamma plays three roles:      
        % gamma is coefficient of relative risk aversion
        % 1/gamma is the intertemporal elasticity of substitution of consumption
        % gamma+1 is the coefficient of relative prudence (Kimball, 1990)

Params.delta=0.03; % rate at which agents discount future
Params.beta=1/(1+Params.delta);
Params.r=0.03; % interest rate on assets

% Table A.1
Params.dj=[0.00058, 0.00061, 0.00062, 0.00064, 0.00065, 0.00067, 0.00069, 0.00070, 0.00072, 0.00075, 0.00078, 0.00082, 0.00086, 0.00091, 0.00098, 0.00105, 0.00115, 0.00128, 0.00144, 0.00161, 0.00180, 0.00200,...
    0.00221, 0.00242, 0.00266, 0.00292, 0.00320, 0.00349, 0.00380, 0.00413, 0.00450, 0.00490, 0.00533,...
    0.00533, 0.00581, 0.00632, 0.00689, 0.00749, 0.00811, 0.00878, 0.00952, 0.01033, 0.01121, 0.01223,...
    0.01332, 0.01455, 0.01590, 0.01730, 0.01874, 0.02028, 0.02203, 0.02404, 0.02623, 0.02863, 0.03128,...
    0.03432, 0.03778, 0.04166, 0.04597, 0.05078, 0.05615, 0.06214, 0.06885, 0.07631, 0.08455, 0.09352,...
    0.10323, 0.11367, 0.12484, 0.13677, 0.14938, 0.16289, 0.17721, 0.19234, 0.20828, 0.22418, 0.23980, 0.25495, 0.26937, 0.28284]; % conditional probability of death
Params.sj=1-Params.dj; % Conditional survival probabilities. Act as a discount rate.


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
% Compare to Fig A.1(a-c): they won't be exactly the same as drop the year
% fixed effects, but eyeballing suggests they are fine.
plot(21:1:100, Params.DeterministicWj.ft1, 21:1:100, Params.DeterministicWj.ft2, 21:1:100, Params.DeterministicWj.ft3)

% Stochastic Wj (Table A.4)
Params.w_rho=[0.955; 0.946; 0.955];
Params.w_sigmasqepsilon=[0.033; 0.025; 0.016];
Params.w_sigmasqu=Params.w_sigmasqepsilon./(1-Params.w_rho.^2);
Params.w_sigmasqupsilon=[0.040; 0.021; 0.014]; % Estimated from PSID but not used in model.
[z1_grid.ft1, pi_z1.ft1]=TauchenMethod(0,Params.w_sigmasqu(1), Params.w_rho(1), n_z(1), Params.q);
[z1_grid.ft2, pi_z1.ft2]=TauchenMethod(0,Params.w_sigmasqu(2), Params.w_rho(2), n_z(1), Params.q);
[z1_grid.ft3, pi_z1.ft3]=TauchenMethod(0,Params.w_sigmasqu(3), Params.w_rho(3), n_z(1), Params.q);

% Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% "This specification ensures that when we compare the certainty case with
% the earnings uncertainty case, we hold the age-conditional means of earnings constant."

% Table A.5
Params.DeterministicMj_working=[5.373+0.073*WorkingAgeVec-0.753*((WorkingAgeVec.^2)/1000); 4.749+0.106*WorkingAgeVec-1.084*((WorkingAgeVec.^2)/1000); 4.816+0.109*WorkingAgeVec-1.090*((WorkingAgeVec.^2)/1000)];
Params.DeterministicMj_retired=[14.441-0.200*RetiredAgeVec+1.288*((RetiredAgeVec.^2)/1000); 11.371-0.101*RetiredAgeVec+0.540*((RetiredAgeVec.^2)/1000); 9.553-0.054*RetiredAgeVec-0.297*((RetiredAgeVec.^2)/1000)];
Params.DeterministicMj=[Params.DeterministicMj_working, Params.DeterministicMj_retired]; % is in logs
Params.m_rho=0.901*ones(3,1);
Params.m_sigmasqepsilon=[0.175; 0.156; 0.153];
Params.m_sigmasqomega=[0.220; 0.220; 0.220];
[z2_grid.ft1, pi_z2.ft1]=TauchenMethod(0,Params.m_sigmasqepsilon(1)/(1-Params.m_rho(1)^2), Params.m_rho(1), n_z(2), Params.q);
[z2_grid.ft2, pi_z2.ft2]=TauchenMethod(0,Params.m_sigmasqepsilon(2)/(1-Params.m_rho(2)^2), Params.m_rho(2), n_z(2), Params.q);
[z2_grid.ft3, pi_z2.ft3]=TauchenMethod(0,Params.m_sigmasqepsilon(3)/(1-Params.m_rho(3)^2), Params.m_rho(3), n_z(2), Params.q);

% Consumption Floor
Params.Cbar=7000; % (middle of pg. 111)

%% Grids
z_grid.ft1=[z1_grid.ft1; z2_grid.ft1];
z_grid.ft2=[z1_grid.ft2; z2_grid.ft2];
z_grid.ft3=[z1_grid.ft3; z2_grid.ft3];

pi_z.ft1=kron(pi_z1.ft1, pi_z2.ft1);
pi_z.ft2=kron(pi_z1.ft2, pi_z2.ft2);
pi_z.ft3=kron(pi_z1.ft3, pi_z2.ft3);

a_grid=linspace(0,maxa,n_a)'; % Could probably do better by adding more grid near zero

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(aprime,a,W_z1,M_z2,gamma,r,Cbar,DeterministicWj, DeterministicMj) HubbardSkinnerZeldes1994_ReturnFn(aprime,a,W_z1,M_z2,gamma,r,Cbar,DeterministicWj, DeterministicMj)
ReturnFnParamNames={'gamma','r','Cbar','DeterministicWj', 'DeterministicMj'}; %It is important that these are in same order as they appear in 'HubbardSkinnerZeldes1994_ReturnFn'

%% Now solve the value function iteration problem

vfoptions.verbose=0;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(0,n_a,n_z,N_j,N_i, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
toc

%% 







