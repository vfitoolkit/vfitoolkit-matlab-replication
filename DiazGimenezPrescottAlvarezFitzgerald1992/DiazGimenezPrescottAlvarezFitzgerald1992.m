% Díaz-Giménez, Prescott, Alvarez & Fitzgerald (1992) - Banking in computable general equilibrium economy

Experiment=0; % 0 for baseline, 1 and 2 for the respective experiment numbers

%% Setup

Params.Experiment=Experiment;

%Real wage, w(s,z)  %calibrated model value is w=[0.1250,0.0400,0,0]';
Params.w1=0.125;
Params.w2=0.04;
Params.w3=0;
Params.w4=0;

%Now we declare all of the other relevant parameters
%Declare and assign values to all the parameters which are calculated
%directly from the calibration (contents of tables on pg 551)
%Calibrate household parameters
%Preferences
Params.alpha=0.3330; %Private consumption share
Params.alpha_k=0.1080; %Capital service share
Params.psi=4.0000; %Risk aversion
Params.beta=0.9994; %Time discount factor
Params.delta_g=0.0104; %Public consumption factor
Params.tau=2.2200; %Productive time
Params.delta_r=0.2100; %Retirees' constant
%Technology
Params.mew=0.00625; %Maintainence cost                     %%%%Note: Value in table in paper is incorrect. Paper reports annual, not eighth-of-year, value.
Params.gamma=0.02625; %Rental service coefficient       %%%% Pg., 550: "we set gamma to twice the sum of the real borrowing rate (i_l-e) and the maintainence cost (mew)"
Params.phi=0.9000; %1-phi = 0.1000; Disinvestment cost
%w the Real wage, w(s,z) is declared above
Params.phi_1=0.9000; %Prob of newborn being: Type 1
Params.phi_2=0.1000; %                       Type 2
%Calibrated bank and government parameters
%Per unit banking costs
Params.eta_d=0.0011875; %Deposits
Params.eta_l=0.0056250; %Loans
%Government policy declared above
Params.rho=0.01;   %reserve requirement(z)
Params.theta=0.2; %Tax rate on labour & interest income(z)
Params.iota=0.00625;    %Nominal interest rate on T-bill(z) 
Params.epsilon_inflation=1.005;    %Inflation rate process
Params.omega=0.02; % Welfare tranfers to Indigent retirees (are zero to everyone else)

%calculate a few of the other parameters directly
Params.i_l=Params.iota+Params.eta_l; %interest rate on loans
Params.i_d=(1-Params.rho)*Params.iota-Params.eta_d; %interest rate on deposits
Params.e=Params.epsilon_inflation-1; %pricing/inflation process on reserves


%% Create grids for s,z,a (A & K), and transition matrices for s & z
n_A=1000; %800 is the original number taken from tech append II, pg 26. Note that n_A also determines the top of the grid.
A_grid=linspace(-2.7,(2.7/245)*(n_A-1)-2.7,n_A)'; %so we n_A points
[~,n_A_zero]=min(abs(A_grid-0)); % find index that corresponds to zero assets (the minus zero is redundant here, but used to make exposition clearer)
A_grid(n_A_zero)=0; % With the original 800 grid points this line is redundant as the grid was designed to include 0, but to use larger grid size we have to ensure that there is still a zero value
K_grid=[0;3]; %household capital stocks (pg 550) (house + consumer durables + small business capital)
n_K=2;

N_grid=[0;1]; %to use for calculating the labour from it's index


%Calibrated household idiosyncratic transtion probablities
pi_s=[0.9593, 0.0369, 0.0038, 0.0000 ; 0.3317, 0.6645, 0.0038, 0.0000 ; 0.0000, 0.0000, 0.9869, 0.0131 ; 0.0000, 0.0000, 0.0000, 1.0000];
s_grid=(1:1:4)';  %individual states: 1&2 working age, 3 retirement, 4 death
% transition matrix for s, the individual specific states row index is this state s, column index is next state s'
%Economy-wide shock transition probabilities
if Experiment==2 %Note that since no results from this appear in any Tables of Figures of the paper it does not actually end up being used as part of this replication.
    pi_z=[0.5,0.5; 0.3,0.7]; %transition matrix for z, row index is z, column index is z'
    z_grid=[1;2];
else
    pi_z=1;
    z_grid=1;
    w_sz=[Params.w1;Params.w2;Params.w3;Params.w4]; % used for welfare evaluations
    sigma_sz=[1;1;1;0]; % used for welfare evaluations
end


%Note that from the point of view of the value function, there is no difference between s & z, so we combine them.
sz_grid=[s_grid;z_grid];
pi_sz=kron(pi_s,pi_z);
n_s=length(s_grid);
n_z=length(z_grid);
n_sz=[n_s,n_z];

a_grid=[A_grid;K_grid];
n_a=[length(A_grid),length(K_grid)];

d_grid=N_grid;
n_d=length(N_grid);

DiscountFactorParamNames={'beta'};

ReturnFn=@(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) DiazGimenezPrescottAlvarezFitzgerald1992_ReturnFn(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment);
ReturnFnParamNames={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};

%% Solve the model
tic;
vfoptions.policy_forceintegertype=1;
[V, Policy]=ValueFnIter_Case1(n_d,n_a,n_sz,d_grid,a_grid,sz_grid, pi_sz, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions); %
time=toc;

fprintf('Time to solve the value function iteration was %8.2f seconds. \n', time)

%% Generate some output following what is reported in Diaz-Gimenez, Prescott, Alvarez & Fitzgerald (1992)

% Not part of replication but let's look at the policy functions
figure(1)
subplot(2,3,1); plot(A_grid,N_grid(Policy(1,:,1,1)),A_grid,N_grid(Policy(1,:,1,2)),A_grid,N_grid(Policy(1,:,1,3)),A_grid,N_grid(Policy(1,:,1,4)))
subplot(2,3,2); plot(A_grid,A_grid(Policy(2,:,1,1)),A_grid,A_grid(Policy(2,:,1,2)),A_grid,A_grid(Policy(2,:,1,3)),A_grid,A_grid(Policy(2,:,1,4)))
subplot(2,3,3); plot(A_grid,K_grid(Policy(3,:,1,1)),A_grid,K_grid(Policy(3,:,1,2)),A_grid,K_grid(Policy(3,:,1,3)),A_grid,K_grid(Policy(3,:,1,4)))

subplot(2,3,4); plot(A_grid,N_grid(Policy(1,:,2,1)),A_grid,N_grid(Policy(1,:,2,2)),A_grid,N_grid(Policy(1,:,2,3)),A_grid,N_grid(Policy(1,:,2,4)))
subplot(2,3,5); plot(A_grid,A_grid(Policy(2,:,2,1)),A_grid,A_grid(Policy(2,:,2,2)),A_grid,A_grid(Policy(2,:,2,3)),A_grid,A_grid(Policy(2,:,2,4)))
subplot(2,3,6); plot(A_grid,K_grid(Policy(3,:,2,1)),A_grid,K_grid(Policy(3,:,2,2)),A_grid,K_grid(Policy(3,:,2,3)),A_grid,K_grid(Policy(3,:,2,4)))


%% First, create the agents stationary distribution
% Requires a customized StationaryDist code as explained below.
StationaryDist=DiazGimenezPrescottAlvarezFitzgerald1992_StationaryDist(Policy,Params,n_d,n_a,n_s,n_z,n_sz,n_A_zero,pi_sz);
% % Simulation of the model is highly non-standard due to the combination of
% % stochastic death in infinitely lived agents together with their not
% % caring about future generations. [Normally would either use
% % finite-lifetime, or if using infinite-lived they would be dynasties
% % that care about future generations.]
% 
% % Note that for n_z=1 this stationary distribution is also the stationary
% % distribution in which all the agents actually find themselves. If n_z>1
% % then it does not represent where all the agents actually at any point in
% % time (which varies with z and must be simulated for a given path of z)
% % but does still represent the asymptotic distribution of agents and so can
% % be used to find asymptotic values of many model statistics of interest.
% 
% % The following makes use of some 'internal functions' of the VFI Toolkit
% % to deal with the non-standard agent distribution. Most of it is simply a
% % minor modification of contents of StationaryDist_Case1().
% simoptions.parallel=2;
% simoptions.tolerance=10^(-9);
% simoptions.maxit=5*10^4;
% PolicyKron=KronPolicyIndexes_Case1(Policy, n_d, n_a, n_sz,simoptions);
% 
% % Create a stationary dist. Am simply setting it up as if all newborns.
% StationaryDist=zeros([n_a,n_sz],'gpuArray');
% StationaryDist(n_A_zero,1,1,:)=Params.phi_1/prod(n_z);
% StationaryDist(n_A_zero,1,2,:)=Params.phi_2/prod(n_z);
% 
% N_a=prod(n_a);
% N_sz=prod(n_sz);
% StationaryDistKron=reshape(StationaryDist,[N_a*N_sz,1]);
% 
% % simoptions.parallel==2 % Using the GPU
% optaprime=reshape(PolicyKron(2,:,:),[1,N_a*N_sz]);
% 
% Ptemp=zeros(N_a,N_a*N_sz,'gpuArray');
% Ptemp(optaprime+N_a*(gpuArray(0:1:N_a*N_sz-1)))=1;
% Ptran=(kron(pi_sz',ones(N_a,N_a,'gpuArray'))).*(kron(ones(N_sz,1,'gpuArray'),Ptemp));
% 
% StationaryDistKronOld=zeros(N_a*N_sz,1,'gpuArray');
% SScurrdist=sum(abs(StationaryDistKron-StationaryDistKronOld));
% SScounter=0;
% 
% while SScurrdist>simoptions.tolerance && (100*SScounter)<simoptions.maxit
%     
%     for jj=1:100
%         StationaryDistKron=Ptran*StationaryDistKron; %No point checking distance every single iteration. Do 100, then check.
%         % NON-STANDARD PART: reallocate the dead to newborns
%         StationaryDistKron=reshape(StationaryDistKron,[N_a,N_sz]); % Could use cleaver indexing, but feeling lazy, so just going to reshape and then un-reshape
%         for z_c=1:n_z
%             MassOfDead=sum(StationaryDistKron(:,n_s*z_c));
%             StationaryDistKron(n_A_zero,1+n_s*(z_c-1))=StationaryDistKron(n_A_zero,1+n_s*(z_c-1))+Params.phi_1*MassOfDead; % want n_A_zero and K==1, but this is anyway just n_A_zero
%             StationaryDistKron(n_A_zero,2+n_s*(z_c-1))=StationaryDistKron(n_A_zero,2+n_s*(z_c-1))+Params.phi_2*MassOfDead;
%             StationaryDistKron(:,n_s*z_c)=0;
%         end
%         StationaryDistKron=reshape(StationaryDistKron,[N_a*N_sz,1]);
%     end
%     
%     StationaryDistKronOld=StationaryDistKron;
%     StationaryDistKron=Ptran*StationaryDistKron;
%     % NON-STANDARD PART: reallocate the dead to newborns
%     StationaryDistKron=reshape(StationaryDistKron,[N_a,N_sz]); % Could use cleaver indexing, but feeling lazy, so just going to reshape and then un-reshape
%     for z_c=1:n_z
%         MassOfDead=sum(StationaryDistKron(:,n_s*z_c));
%         StationaryDistKron(n_A_zero,1+n_s*(z_c-1))=StationaryDistKron(n_A_zero,1+n_s*(z_c-1))+Params.phi_1*MassOfDead; % want n_A_zero and K==1, but this is anyway just n_A_zero
%         StationaryDistKron(n_A_zero,2+n_s*(z_c-1))=StationaryDistKron(n_A_zero,2+n_s*(z_c-1))+Params.phi_2*MassOfDead;
%         StationaryDistKron(:,n_s*z_c)=0;
%     end
%     StationaryDistKron=reshape(StationaryDistKron,[N_a*N_sz,1]);
% 
%     SScurrdist=sum(abs(StationaryDistKron-StationaryDistKronOld));    
%     
%     SScounter=SScounter+1;
%     if rem(SScounter,50)==0
%         SScounter
%         SScurrdist
%     end
% end
% 
% StationaryDist=reshape(StationaryDistKron,[n_a,n_sz]);
% 
% sum(reshape(StationaryDist,[N_a,N_sz]),1)

%% Now we have StationaryDist, time to start replicating output.
save ./SavedOutput/DPAF1992.mat StationaryDist V Policy

%% Replicate Tables 7a & 7b
if Experiment==0
    %Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
    %pi_s (used to calculate the integral across the SS dist fn of whatever
    %functions you define here)
    FnsToEvaluateParamNames=struct();
    FnsToEvaluateParamNames(1).Names={};
    FnsToEvaluateFn_1 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val) a2_val; % Capital (Housing Stock)
    FnsToEvaluateParamNames(2).Names={'e'};
    FnsToEvaluateFn_2 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,e) (1+e)*a1_val*(a1_val>0); % Deposits
    FnsToEvaluateParamNames(3).Names={'e'};
    FnsToEvaluateFn_3 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,e) -(1+e)*a1_val*(a1_val<0); % Loans
    FnsToEvaluateParamNames(4).Names={'e'};
    FnsToEvaluateFn_4 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,e) (1+e)*a1_val+a2_val; % Net Worth: here just as a double-check as should equal Capital+Deposits-Loans
    % Reserves: equals reserve requirement rho, multiplied by Deposits
    % Gov Debt: equals Deposits-Loans (by banks balance sheet)
    FnsToEvaluate={FnsToEvaluateFn_1, FnsToEvaluateFn_2, FnsToEvaluateFn_3, FnsToEvaluateFn_4};
    AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_sz, d_grid, a_grid,sz_grid, 2);
    
    AggVars
    
    FnsToEvaluateParamNames=struct();
    FnsToEvaluateParamNames(1).Names={'phi'};
    FnsToEvaluateFn_1 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi) a2prime_val*(a2prime_val-a2_val>0)-phi*a2_val*(a2prime_val-a2_val<0); % Investment (Housing investment)
    FnsToEvaluateParamNames(2).Names={'gamma'};
    FnsToEvaluateFn_2 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,gamma) gamma*(a2_val==0); % Imputed rents of non-homeowners
    FnsToEvaluateParamNames(3).Names={'eta_d','eta_l'};
    FnsToEvaluateFn_3 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,eta_d,eta_l) eta_d*a1_val*(a1_val>0)-eta_l*a1_val*(a1_val<0); % Banking Value Added
    FnsToEvaluateParamNames(4).Names={'w1','w2'};
    FnsToEvaluateFn_4 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,w1,w2) w1*(s_val==1)*d_val+w2*(s_val==2)*d_val; % Wage income
    FnsToEvaluateParamNames(5).Names={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};
    FnsToEvaluateFn_5 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) DiazGimenezPrescottAlvarezFitzgerald1992_ConsumptionFn(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment); % Goods Consumption of Homeowners
    FnsToEvaluateParamNames(6).Names={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};
    FnsToEvaluateFn_6 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) (s_val==1||s_val==2)*DiazGimenezPrescottAlvarezFitzgerald1992_ConsumptionFn(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment); % Goods Consumption of workers
    FnsToEvaluateParamNames(7).Names={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};
    FnsToEvaluateFn_7 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) (s_val==1||s_val==2)*(a2_val>0)*DiazGimenezPrescottAlvarezFitzgerald1992_ConsumptionFn(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment); % Goods Consumption of Homeowners
    FnsToEvaluateParamNames(8).Names={};
    FnsToEvaluateFn_8 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val) a2_val*(s_val==1||s_val==2); % Capital (Housing Stock) owned by workers
    FnsToEvaluateParamNames(9).Names={'omega'};
    FnsToEvaluateFn_9 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,omega) omega*(s_val==3); % Pension income
    FnsToEvaluate={FnsToEvaluateFn_1, FnsToEvaluateFn_2, FnsToEvaluateFn_3, FnsToEvaluateFn_4, FnsToEvaluateFn_5, FnsToEvaluateFn_6, FnsToEvaluateFn_7, FnsToEvaluateFn_8, FnsToEvaluateFn_9};
    AggVars2=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_sz, d_grid, a_grid,sz_grid, 2);
    
    AggVars2 %*(2.7/2.24)
    
    % Note: model period is one-eigth of a year and the entries in Table 7b
    % all relate to annual GDP (hence the multiplying by eight in many of
    % the below).
%     GNP=; % should be about 1/1.2054
% GNP = 8*SSvalues_AggVars2(4)
    GNP=2.24/2.7;
    % = 0.8296 % Just force it to take known value for now
    % Wage income: SSvalues_AggVars2(4)*8 = 0.6987
    % Pension income: SSvalues_AggVars2(9)*8 = 0.036
    % Income of banks: 8*Params.eta_l*SSvalues_AggVars(3)+8*Params.eta_d*SSvalues_AggVars(2) = 0.025
    
%     HousingValueAdded= % 15.74 
    BankingValueAdded=8*Params.eta_l*AggVars(3)/GNP+8*Params.eta_d*AggVars(2)/GNP; % 3.01
%     GoodsProducingValueAdded= % 81.25

        Goods=8*AggVars2(5)/GNP; % 50.3
        % Housing=(alpha_k/alpha)*SSvalues_AggVars2(7)*8/GNP; 
%         Housing=8*(Params.gamma-Params.i_l)*SSvalues_AggVars(1)/GNP; % 29.3
            Maintainence=8*Params.mew*AggVars(1)/GNP; % 13.6
%         BankingServices= % 3.1
% 8*(Params.i_l)*SSvalues_AggVars(3)/GNP+8*(Params.i_d)*SSvalues_AggVars(2)/GNP-8*(0.05/8)*(SSvalues_AggVars(2)-SSvalues_AggVars(3))/GNP
% 8*(Params.i_l)*(SSvalues_AggVars(3)/(1-Params.rho))/GNP+8*(Params.i_d)*SSvalues_AggVars(2)/GNP
% 8*(Params.eta_l)*(SSvalues_AggVars(3)/(1-Params.rho))/GNP+8*(Params.eta_d)*SSvalues_AggVars(2)/GNP
% 
% 0.05*(SSvalues_AggVars(2)-SSvalues_AggVars(3))/GNP
%     Consumption=Goods+Housing+BankingServices; % 82.6
    GovPurchases=8*(AggVars2(4)-(AggVars2(5)+AggVars2(1)+AggVars2(3)+Params.mew*AggVars(1)))/GNP; % 16.6 [=nw-(c+(xd-xs)+(detad+letal)+mewk')]
    Investment=8*AggVars2(1)/GNP; %=x_d+x_s % 0.84

%     (alpha_k/alpha)*8*(SSvalues_AggVars2(5)-SSvalues_AggVars2(6))/GNP % Cons by non-owning workers
    
    % Imputed rents from all houses: 8*(Params.gamma-Params.i_l)*SSvalues_AggVars(1)/GNP  =0.3113
    %                                8*Params.gamma*SSvalues_AggVars(1)/GNP = 0.5685
    % Imputed rents from all houses owned by workers: 8*(Params.gamma-Params.i_l)*SSvalues_AggVars2(8)/GNP   =0.3101
    % Housing services consumption of renters: 
    
    FilenameString=['./SavedOutput/LatexInputs/DiazGimenezPrescottAlvarezFitzgerald1992_Tables7a7b.tex'];
    FID = fopen(FilenameString, 'w');
    fprintf(FID, 'Calibrated economys steady-state balance sheet data \\\\ \n');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lr} \\hline \n');
    fprintf(FID, '  & Stock/GNP \\\\ \\hline \n');
    fprintf(FID, ' \\textit{Household sector} & \\\\ \n');
    fprintf(FID, ' \\quad Tangible Capital & %1.2f \\\\ \n', AggVars(1)/GNP);
    fprintf(FID, ' \\quad Deposits & %1.2f \\\\ \n', AggVars(2)/GNP);
    fprintf(FID, ' \\quad Loans & %1.2f \\\\ \n', AggVars(3)/GNP);
    fprintf(FID, ' \\quad Net Worth & %1.2f \\\\ \n', AggVars(4)/GNP);
    fprintf(FID, ' \\textit{Government Sector} & \\\\ \n');
    fprintf(FID, ' \\quad Reserves & %1.2f \\\\ \n', Params.rho*AggVars(2)/GNP);
    fprintf(FID, ' \\quad Debt & %1.2f \\\\ \n', (AggVars(2)-AggVars(3))/GNP);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\vspace{5mm} \n');
    fprintf(FID, 'Calibrated economys steady-state NIPA data \\\\ \n');
    fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lr} \\hline \n');
    fprintf(FID, '  & Percent of GNP \\\\ \\hline \n');
    fprintf(FID, ' \\textit{Value added by Sectors} & \\\\ \n');
    fprintf(FID, ' \\quad Housing &  \\\\ \n'); %MISSING 100*HousingValueAdded
    fprintf(FID, ' \\quad Banking & %1.2f \\\\ \n', 100*BankingValueAdded);
    fprintf(FID, ' \\quad Goods producing & \\\\ \n'); %MISSING 100*GoodsProducingValueAdded
    fprintf(FID, ' \\textit{Production} & \\\\ \n');
    fprintf(FID, ' \\quad Consumption &  \\\\ \n'); % MISSING 100*Consumption
    fprintf(FID, ' \\quad \\quad Goods & %1.2f \\\\ \n', 100*Goods);
    fprintf(FID, ' \\quad \\quad Housing & \\\\ \n'); % MISSING 100*Housing
    fprintf(FID, ' \\quad \\quad \\quad Maintainence & %1.2f \\\\ \n', 100*Maintainence);
    fprintf(FID, ' \\quad \\quad Banking Services &  \\\\ \n'); % MISSING , 100*BankingServices
    fprintf(FID, ' \\quad Government Purchases & %1.2f \\\\ \n', 100*GovPurchases);
    fprintf(FID, ' \\quad Investment & %1.2f \\\\ \n', 100*Investment);
    fprintf(FID, '\\hline \n \\end{tabular*} \n');
    fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
    fprintf(FID, 'Replication of Tables 7a \\& 7b of Diaz-Gimenez, Prescott, Alvarez \\& Fitzgerald (1992) using grid sizes $n_a=%d $, $ n_s=%d $, $ n_z=%d $ (Capital and Labour supply are discrete variables with two possible values each) \\\\ \n', n_A, n_s, n_z);
    fprintf(FID, 'Missing values for Housing value added, Goods producing value added, Consumption, Housing, and Banking Services reflect numbers that could not be replicated. Since original paper does not report explicit formulae for their calculation it is unclear if this was simply an inability to find correct formula or actual error in original paper.');
    fprintf(FID, '}} \\end{minipage}');
    fclose(FID);
end


%% Welfare Comparisons
% Perform some calculations required as part of the Welfare comparison formulae
%[ValueFnofPublicandPrivateConsumption,v1,v2,W]=DiazGimenezPrescottAlvarezFitzgerald1992_Welfare(Params,V,n_A,n_K,n_s,n_z,n_d,n_a,n_sz,sigma_sz,w_sz,pi_sz,A_grid,K_grid,d_grid, a_grid,sz_grid,StationaryDist, Policy);
% 
% % h(s,z): lifetime human capital wealth
% hdist=Inf;
% h=zeros(n_s,n_z);
% while hdist>10^(-9)
%     hold=h;
%     h=sigma_sz.*w_sz+Params.beta*sigma_sz.*(pi_sz*hold); % might need to be pi_sz'
% 
%     hdist=sum(abs(h-hold));    
% end
% 
% % Total Wealth of Household
% W=zeros(n_A,n_K,n_s,n_z);
% for i1=1:n_A
%     for i2=1:n_K
%         for i3=1:n_s
%             for i4=1:n_z
%                 W(i1,i2,i3,i4)=A_grid(i1)+K_grid(i2)+h(i3,i4);
%             end
%         end
%     end
% end
% 
% % v1 is the value function computed earlier
% 
% % v2: value of public consumption (note: this code only works for n_z=1)
% FnsToEvaluateParamNames=struct();
% FnsToEvaluateParamNames(1).Names={'phi'};
% FnsToEvaluateParamNames(2).Names={'eta_d','eta_l'};
% FnsToEvaluateParamNames(3).Names={'w1','w2'};
% FnsToEvaluateParamNames(4).Names={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};
% FnsToEvaluateParamNames(5).Names={};
% FnsToEvaluateFn_1 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi) a2prime_val*(a2prime_val-a2_val>0)-phi*a2_val*(a2prime_val-a2_val<0); % Investment (Housing investment)
% FnsToEvaluateFn_2 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,eta_d,eta_l) eta_d*a1_val*(a1_val>0)-eta_l*a1_val*(a1_val<0); % Banking Value Added
% FnsToEvaluateFn_3 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,w1,w2) w1*(s_val==1)*d_val+w2*(s_val==2)*d_val; % Wage income
% FnsToEvaluateFn_4 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) DiazGimenezPrescottAlvarezFitzgerald1992_ConsumptionFn(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment); % Goods Consumption of Homeowners
% FnsToEvaluateFn_5 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val) a2_val; % Capital (Housing Stock)
% FnsToEvaluate={FnsToEvaluateFn_1, FnsToEvaluateFn_2, FnsToEvaluateFn_3, FnsToEvaluateFn_4, FnsToEvaluateFn_5};
% AggVars3=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_sz, d_grid, a_grid,sz_grid, 2);
% g=8*(AggVars3(3)-(AggVars3(4)+AggVars3(1)+AggVars3(2)+Params.mew*AggVars3(5)));
% 
% v2dist=Inf;
% v2=zeros(n_s,n_z);
% while v2dist>10^(-9)
%     v2old=v2;
%     v2=sigma_sz*Params.delta_g*(g^(Params.alpha*(1-Params.psi)))/(1-Params.psi)+Params.beta*sigma_sz.*(pi_sz*v2old);
%     
%     v2dist=sum(abs(v2-v2old));
% end
% v2=real(v2); % for some reason matlab is making them complex numbers even though imaginary component is 0
% v2=gather(v2);
% 
% % V= v1+v2
% v1=gather(V);
% ValueFnofPublicandPrivateConsumption=zeros(n_A,n_K,n_s,n_z);
% for i1=1:n_A
%     for i2=1:n_K
%         for i3=1:n_s
%             for i4=1:n_z
%                 ValueFnofPublicandPrivateConsumption(i1,i2,i3,i4)=v1(i1,i2,i3,i4)+v2(i3,i4);
%             end
%         end
%     end
% end

%% Table 8
% Uses a substantially simplified and reparametrized version of the model. Hence have
% just put it in an seperate command to be called.
Table8=DiazGimenezPrescottAlvarezFitzgerald1992_Table8;

FilenameString=['./SavedOutput/LatexInputs/DiazGimenezPrescottAlvarezFitzgerald1992_Table8.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, 'Benefits of Switching to a policy of a zero after-tax real return on deposits \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \\hline \n');
fprintf(FID, 'Current after-tax  & \\multicolumn{6}{c}{Benefits as percent of wealth} \\\\ \\cline{2-7} \n');
fprintf(FID, 'real return        & \\multicolumn{3}{c}{Measured at Mean} & \\multicolumn{3}{c}{Measured across Stationary Dist} \\\\ \\hline \n');
fprintf(FID, 'on deposits        & $S$ & $S_1$ & $S_2$ & $M$ & $M_1$ & $M_2$ \\\\ \\hline \n');
fprintf(FID, ' %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', Table8(1,:));
fprintf(FID, ' %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', Table8(2,:));
fprintf(FID, ' %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', Table8(3,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Replication of Table 8 of Diaz-Gimenez, Prescott, Alvarez \\& Fitzgerald (1992) using grid sizes $n_a=%d $, $ n_s=2 $ \\\\ \n', n_A);
fprintf(FID, 'Original paper incorrectly describes the first measure as evaluated at the steady-state, it is in fact being evaluated at the mean of the stationary distribution. ');
fprintf(FID, 'Note that Table 8 produces results relating not directly to the model of the paper itself, but to the model of \\cite*{ImrohorogluPrescott1991}, which is a subcase of the model, hence the different grid sizes used for this Table compared to others from this paper.');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Table 9
Table9=zeros(3,7);
inflation_vec=[0.04,0.05,0.06]./8; % 4,5 & 6 percent inflation for Table 9. Model period is one-eigth of year.

% Generate objects needed to calculate rows of Table 8 (the reference case and the three rows)
for ii=1:3 
    Params.epsilon_inflation=1+inflation_vec(ii);    %Inflation rate process
    Params.e=Params.epsilon_inflation-1; %pricing/inflation process on reserves
    
    V0=ones([n_a,n_sz]);
    vfoptions.policy_forceintegertype=1;
    [V, Policy]=ValueFnIter_Case1(n_d,n_a,n_sz,d_grid,a_grid,sz_grid, pi_sz, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions); %
    
    StationaryDist=DiazGimenezPrescottAlvarezFitzgerald1992_StationaryDist(Policy,Params,n_d,n_a,n_s,n_z,n_sz,n_A_zero,pi_sz);
    
    [ValueFnofPublicandPrivateConsumption,v1,v2,W]=DiazGimenezPrescottAlvarezFitzgerald1992_Welfare(Params,V,n_A,n_K,n_s,n_z,n_d,n_a,n_sz,sigma_sz,w_sz,pi_sz,A_grid,K_grid,d_grid, a_grid,sz_grid,StationaryDist, Policy);

    Outputs(ii).V=gather(V);
    Outputs(ii).Policy=gather(Policy);
    Outputs(ii).StationaryDist=gather(StationaryDist);
    Outputs(ii).ValueFnofPublicandPrivateConsumption=gather(ValueFnofPublicandPrivateConsumption);
    Outputs(ii).v1=gather(v1);
    Outputs(ii).v2=gather(v2);
    Outputs(ii).W=gather(W);
    
end

% Now calculate the actual entries
for ii=2:3
    Lambda=(Outputs(1).ValueFnofPublicandPrivateConsumption ./ Outputs(ii).ValueFnofPublicandPrivateConsumption).^(1/(Params.alpha*(1-Params.psi)));
    
    % Some elements of Lambda have ended up 'nan'. Checked and they all
    % correspond to -Inf/-Inf. Replace them with 1 (so Lambda-1 will equal
    % zero and they will have no impact on M)
    Lambda(isnan(Lambda))=1;
    
    % n_z=1, so can ignore this aspect of formula for M (pg 554) (y in
    % formula on pg 554 is the stationary distribution for the pi0 policy)
    M=sum(sum(sum(Outputs(ii).W.*(Lambda-1).*Outputs(ii).StationaryDist)))/sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));
    
    M1factor=(Outputs(1).v1-Outputs(ii).v1)./(Outputs(1).ValueFnofPublicandPrivateConsumption-Outputs(ii).ValueFnofPublicandPrivateConsumption);
    M2factor=(ones(n_A,n_K).*shiftdim(Outputs(1).v2-Outputs(ii).v2,-2))./(Outputs(1).ValueFnofPublicandPrivateConsumption-Outputs(ii).ValueFnofPublicandPrivateConsumption);
    
    % Same problem of dividing by -Inf that had with Lambda in M occours in
    % M1factor and M2factor. This time replace them by zeros.
    M1factor(isnan(M1factor))=0;
    M2factor(isnan(M2factor))=0;
    
    M1=sum(sum(sum(Outputs(ii).W.*M1factor.*(Lambda-1).*Outputs(ii).StationaryDist)))/sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));
    M2=sum(sum(sum(Outputs(ii).W.*M2factor.*(Lambda-1).*Outputs(ii).StationaryDist)))/sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));

    
    Table9(ii,5)=M;
    Table9(ii,6)=M1;
    Table9(ii,7)=M2;
end

Table9(:,1)=inflation_vec'.*8;
Table9(:,2)=((1-Params.theta)*Params.i_d-inflation_vec')*8;
Table9(:,3)=ones(3,1)*inflation_vec(1).*8;
Table9(:,4)=((1-Params.theta)*Params.i_d-ones(3,1)*inflation_vec(1))*8;

Table9=Table9*100;

FilenameString=['./SavedOutput/LatexInputs/DiazGimenezPrescottAlvarezFitzgerald1992_Table9.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, 'Welfare benefits of switching to a policy of less-negative after-tax real return on deposits \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \\hline \n');
fprintf(FID, ' & &  & & \\multicolumn{3}{c}{Benefits} \\\\ \n');
fprintf(FID, ' \\multicolumn{2}{c}{Current policy} & \\multicolumn{2}{c}{New policy} & \\multicolumn{3}{c}{(percent of wealth)} \\\\ \\cline{1-2} \\cline{3-4} \\cline{5-7} \n');
fprintf(FID, 'Inflation & After-tax & Inflation & After-tax & & &  \\\\ \n');
fprintf(FID, 'on deposits & real return & rate & real return & Total & Private & Public \\\\ \\hline \n');
fprintf(FID, ' %d\\%% & %1.2f\\%% & %d\\%% & %1.2f\\%% &  &  &  \\\\ \n', Table9(1,1:4));
fprintf(FID, ' %d & %1.2f & %d & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', Table9(2,:));
fprintf(FID, ' %d & %1.2f & %d & %1.2f & %1.2f & %1.2f & %1.2f \\\\ \n', Table9(3,:));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Replication of Table 9 of Diaz-Gimenez, Prescott, Alvarez \\& Fitzgerald (1992) using grid sizes $n_a=%d $, $ n_s=%d $. \n', n_A, n_s);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

