% Díaz-Giménez, Prescott, Alvarez & Fitzgerald (1992) - Banking in
% computable general equilibrium economy

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
Params.mew=0.00625; %Maintainence cost                     %%%%Note: Value in table in paper is incorrect.
Params.gamma=0.02625; %Rental service coefficient
Params.phi=0.9000; %1-phi_lc = 0.1000; Disinvestment cost
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
n_A=800; %number taken from tech append II, pg 26
A_grid=linspace(-2.7,(2.7/245)*(n_A-1)-2.7,n_A)'; %so we n_A points
[~,n_A_zero]=min(abs(A_grid-0)); % find index that corresponds to zero assets (the minus zero is redundant here, but used to make exposition clearer)
K_grid=[0;3]; %household capital stocks (pg 550) (house + consumer durables + small business capital)

N_grid=[0;1]; %to use for calculating the labour from it's index


%Calibrated household idiosyncratic transtion probablities
pi_s=[0.9593, 0.0369, 0.0038, 0.0000 ; 0.3317, 0.6645, 0.0038, 0.0000 ; 0.0000, 0.0000, 0.9869, 0.0131 ; 0.0000, 0.0000, 0.0000, 1.0000];
s_grid=(1:1:4)';  %individual states: 1&2 working age, 3 retirement, 4 death
% transition matrix for s, the individual specific states row index is this state s, column index is next state s'
%Economy-wide shock transition probabilities
if Experiment==2
    pi_z=[0.5,0.5; 0.3,0.7]; %transition matrix for z, row index is z, column index is z'
    z_grid=[1;2];
else
    pi_z=1;
    z_grid=1;
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

ReturnFn=@(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) DiazGimenezPrescottAlvarezFitzgerald1992_ReturnFn(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment)
ReturnFnParamNames={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};

%% Solve the model
tic;
V0=ones([n_a,n_sz]);
[V, Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_sz,d_grid,a_grid,sz_grid, pi_sz, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);
time=toc;

fprintf('Time to solve the value function iteration was %8.2f seconds. \n', time)


%% Generate some output following what is reported in Diaz-Gimenez, Prescott, Alvarez & Fitzgerald (1992)
% NOT YET DONE


% Simulation of the model is highly non-standard due to the combination of
% stochastic death in infinitely lived agents together with their not
% caring about future generations. [Normally would either use
% finite-lifetime, or if using infinite-lived they would be dynasties
% that care about future generations.]

% The following makes use of some 'internal functions' of the VFI Toolkit
% to deal with the non-standard agent distribution. Most of it is simply a
% minor modification of contents of StationaryDist_Case1().
simoptions.tolerance=10^(-9);
simoptions.maxit=5*10^4;
PolicyKron=KronPolicyIndexes_Case1(Policy, n_d, n_a, n_sz,simoptions);

% Create a stationary dist. Am simply setting it up as if all newborns.
StationaryDist=zeros([n_a,n_sz],'gpuArray');
StationaryDist(n_A_zero,1,1,:)=Params.phi_1/prod(n_z);
StationaryDist(n_A_zero,1,2,:)=Params.phi_2/prod(n_z);

N_a=prod(n_a);
N_sz=prod(n_sz);
StationaryDistKron=reshape(StationaryDist,[N_a*N_sz,1]);

% simoptions.parallel==2 % Using the GPU
optaprime=reshape(PolicyKron(2,:,:),[1,N_a*N_sz]);

Ptemp=zeros(N_a,N_a*N_sz,'gpuArray');
Ptemp(optaprime+N_a*(gpuArray(0:1:N_a*N_sz-1)))=1;
Ptran=(kron(pi_sz',ones(N_a,N_a,'gpuArray'))).*(kron(ones(N_sz,1,'gpuArray'),Ptemp));

StationaryDistKronOld=zeros(N_a*N_sz,1,'gpuArray');
SScurrdist=sum(abs(StationaryDistKron-StationaryDistKronOld));
SScounter=0;

while SScurrdist>simoptions.tolerance && (100*SScounter)<simoptions.maxit
    
    for jj=1:100
        StationaryDistKron=Ptran*StationaryDistKron; %No point checking distance every single iteration. Do 100, then check.
        % NON-STANDARD PART: reallocate the dead to newborns
        StationaryDistKron=reshape(StationaryDistKron,[N_a,N_sz]); % Could use cleaver indexing, but feeling lazy, so just going to reshape and then un-reshape
        for z_c=1:n_z
            MassOfDead=sum(StationaryDistKron(:,n_s*z_c));
            StationaryDistKron(n_A_zero,1+n_s*(z_c-1))=StationaryDistKron(n_A_zero,1+n_s*(z_c-1))+Params.phi_1*MassOfDead; % want n_A_zero and K==1, but this is anyway just n_A_zero
            StationaryDistKron(n_A_zero,2+n_s*(z_c-1))=StationaryDistKron(n_A_zero,2+n_s*(z_c-1))+Params.phi_2*MassOfDead;
            StationaryDistKron(:,n_s*z_c)=0;
        end
        StationaryDistKron=reshape(StationaryDistKron,[N_a*N_sz,1]);
    end
    
    StationaryDistKronOld=StationaryDistKron;
    StationaryDistKron=Ptran*StationaryDistKron;
    % NON-STANDARD PART: reallocate the dead to newborns
    StationaryDistKron=reshape(StationaryDistKron,[N_a,N_sz]); % Could use cleaver indexing, but feeling lazy, so just going to reshape and then un-reshape
    for z_c=1:n_z
        MassOfDead=sum(StationaryDistKron(:,n_s*z_c));
        StationaryDistKron(n_A_zero,1+n_s*(z_c-1))=StationaryDistKron(n_A_zero,1+n_s*(z_c-1))+Params.phi_1*MassOfDead; % want n_A_zero and K==1, but this is anyway just n_A_zero
        StationaryDistKron(n_A_zero,2+n_s*(z_c-1))=StationaryDistKron(n_A_zero,2+n_s*(z_c-1))+Params.phi_2*MassOfDead;
        StationaryDistKron(:,n_s*z_c)=0;
    end
    StationaryDistKron=reshape(StationaryDistKron,[N_a*N_sz,1]);

    SScurrdist=sum(abs(StationaryDistKron-StationaryDistKronOld));    
    
    SScounter=SScounter+1;
    if rem(SScounter,50)==0
        SScounter
        SScurrdist
    end
end

StationaryDist=reshape(StationaryDistKron,[n_a,n_sz]);

%% Now we have StationaryDist, time to start replicating output.

%% Replicate Tables 7a & 7b
if Experiment==0
    
end










