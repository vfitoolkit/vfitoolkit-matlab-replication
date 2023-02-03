% Kulish, Kent & Smith (2010) - Aging, Retirement and Savings: A General Equilibrium Analysis
%
% Note: the baseline model is only 60 periods, but an extension has 70
% periods. Here we model 70 periods, and just make the last 10 irrelevant
% in the baseline setup.
%
% Note: There is NO endogenous retirement "decision" in this model.
% There is endogenous labor supply and at later ages the agents will choose
% to supply zero labor (or more precisely choose leisure=1, which is the
% actual decision in the model) because of a changing preference for leisure, denoted v(j,J). 
% But there is no explicit decision to retire.

n_d=51;
n_assets=701;
n_z=0; % No idiosyncratic shocks

N_j=70; % Needs to be up to 70, but in baseline it is effectively 60

CreateFigures=0

% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

%% Parameters

Params.J=60; % Number of years of life (at most N_j; KKS2010 call this T)
Params.agej=1:1:N_j;
% To enforce that agents live only up to age 60 I use a discount factor so
% they don't care after this
Params.survival=[ones(1,Params.J),zeros(1,N_j-Params.J)];

Params.n=0.012; % Growth rate of cohort size

% Preferences
Params.beta=0.97; % Discount factor
Params.rho=3.5;

% Deterministic earnings profile (KKS2010 call this es)
Params.ej=2.5*[1.0000000000,1.0375471841, 1.0751431228, 1.1126918448, 1.1500940792, 1.1872476552, 1.2240479353, 1.2603882795, 1.2961605375, 1.3312555663, 1.3655637700, 1.3989756576, 1.4313824160, 1.4626764928, 1.4927521862, 1.5215062353, 1.5488384090, 1.5746520864, 1.5988548252, 1.6213589134, 1.6420818991, 1.6609470952, 1.6778840524, 1.6928289993, 1.7057252429, 1.7165235285, 1.7251823538, 1.7316682356, 1.7359559260, 1.7380285767, 1.7378778496, 1.7355039725, 1.7309157401, 1.7241304590, 1.7151738387, 1.7040798284, 1.6908904014, 1.6756552896, 1.6584316695, 1.6392838031, 1.6182826364, 1.5955053603, 1.5710349353, 1.5449595871, 1.5173722747, 1.4883701373, 1.4580539237, 1.4265274094, 1.3938968055, 1.3602701653, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914, 1.3257567914]; 
% The vector comes from "Human capital profile.xls" which Mariano Kulish sent me by email. 
% Because Figure 2 has first value of 2.5, I assume this should be multiplied by 2.5 (result 
% looks just like Figure 2, so fairly confidend of this).
% Alternative
Params.ej_alternative=[Params.ej(1:50),linspace(Params.ej(end),2,70-50)]; % I just eyeballed the tail of ej_alternative, but it is anyway not actually used for anything in the paper (except Figure 2)

Params.g_e=0; % Deterministic growth rate of labor efficiency units (needed for Section 4.6; KKS2010 call this n_e, but that would confuse it with VFI Toolkit notation)

% v(j,J) determines how utility of leisure changes with age and is
% controlled by three parameters
Params.b1=3;
Params.b2=0.7;
Params.b3=0.03;

% Production function
Params.A=0.4;
Params.alpha=0.45;
Params.delta=0.052;

% Initial guess for the general eqm parameter
Params.r=0.03;


%% Create Table 1 and Figure 2, the calibrated parameters
FID = fopen('./SavedOutput/LatexInputs/KulishKentSmith2010_Table1.tex', 'w');
fprintf(FID, '\\begin{tabular}{lll} \n \\hline \n');
fprintf(FID, 'Variable & Description & Value \\\\ \\hline \n');
fprintf(FID, 'T (years) & Life expectancy & %i \\\\ \n',Params.J);
fprintf(FID, 'n & Cohort growth rate & %8.3f \\\\ \n',Params.n);
fprintf(FID, '$e_j$ & Human capital profile & See Figure 1 \\\\ \n');
fprintf(FID, '$b_1$ & Parameter of $v(j,T)$ & %8.1f \\\\ \n',Params.b1);
fprintf(FID, '$b_2$ & Parameter of $v(j,T)$ & %8.1f \\\\ \n',Params.b2);
fprintf(FID, '$b_3$ & Parameter of $v(j,T)$ & %8.2f \\\\ \n',Params.b3);
fprintf(FID, '$\\beta$ & Human capital profile & %8.2f \\\\ \n',Params.beta);
fprintf(FID, '$\\rho$ & Human capital profile & %8.1f \\\\ \n',Params.rho);
fprintf(FID, 'A & Human capital profile & %8.1f \\\\ \n',Params.A);
fprintf(FID, '$1-\\alpha$ & Human capital profile & %8.2f \\\\ \n',1-Params.alpha);
fprintf(FID, '$\\delta$ & Human capital profile & %8.3f \\\\ \n',Params.delta);
fprintf(FID, '\\hline \n \\end{tabular} \n');
fclose(FID);

if CreateFigures==1
    fig1=figure(1);
    plot(1:1:N_j,Params.ej_alternative,1:1:N_j,Params.ej)
    legend('Alternative','Benchmark')
    xlabel('Age-years')
    title('Human Capital Profile')
    saveas(fig1,'./SavedOutput/Graphs/KulishKentSmith2010_Figure2.png')
end

%% Check that v(j,J) is set up correctly by plotting Figure 1
% Note: have had to just guess what b2star and b3star should be

if CreateFigures==1
    fig2=figure(2);
    
    % This graph does not use the actual model parameterization
    Params.b1graph=1.25;
    Params.b2graph=0.7;
    Params.b3graph=0.05;
    Params.Jgraph=60;
    
    % Create the standard for the graphs
    v_jJ=ones(1,N_j);
    for jj=1:Params.Jgraph
        v_jJ(jj)=(Params.b1graph*jj/Params.Jgraph)*normcdf(jj,Params.b2graph*Params.Jgraph,Params.b3graph*Params.Jgraph);
    end
    
    Params.b1star=1; % This can be visually inferred from the graph of KKS2010
    v_jJ1=ones(1,N_j);
    for jj=1:Params.Jgraph
        v_jJ1(jj)=(Params.b1star*jj/Params.Jgraph)*normcdf(jj,Params.b2graph*Params.Jgraph,Params.b3graph*Params.Jgraph);
    end
    subplot(2,2,1); plot(Params.agej(1:60),v_jJ(1:60),Params.agej(1:60),v_jJ1(1:60))
    xlim([1,70])
    
    Params.b2star=0.53; % This is just a guess (will move the mean to roughly age 32)
    v_jJ2=ones(1,N_j);
    for jj=1:Params.Jgraph
        v_jJ2(jj)=(Params.b1graph*jj/Params.Jgraph)*normcdf(jj,Params.b2star*Params.Jgraph,Params.b3graph*Params.Jgraph);
    end
    subplot(2,2,2); plot(Params.agej(1:60),v_jJ(1:60),Params.agej(1:60),v_jJ2(1:60))
    xlim([1,70])
    
    Params.b3star=0.1; % This is just a guess
    v_jJ3=ones(1,N_j);
    for jj=1:Params.Jgraph
        v_jJ3(jj)=(Params.b1graph*jj/Params.Jgraph)*normcdf(jj,Params.b2graph*Params.Jgraph,Params.b3star*Params.Jgraph);
    end
    subplot(2,2,3); plot(Params.agej(1:60),v_jJ(1:60),Params.agej(1:60),v_jJ3(1:60))
    xlim([1,70])
    
    Params.Jstar=70;
    v_jJ4=ones(1,N_j);
    for jj=1:Params.Jstar
        v_jJ4(jj)=(Params.b1graph*jj/Params.Jstar)*normcdf(jj,Params.b2graph*Params.Jstar,Params.b3graph*Params.Jstar);
    end
    subplot(2,2,4); plot(Params.agej(1:60),v_jJ(1:60),Params.agej(1:70),v_jJ4(1:70))
    xlim([1,70])
    
    saveas(fig2,'./SavedOutput/Graphs/KulishKentSmith2010_Figure1.png')
end


%%
d_grid=linspace(0,1,n_d)'; % Grid for leisure

assets_grid=50*(linspace(0,1,n_assets).^3)'; % Grid for assets

retirement_grid=[0,1]; % 1=retired

n_a=n_assets;
a_grid=assets_grid;

z_grid=[];
pi_z=1;


%% Now, create the return function
DiscountFactorParamNames={'beta','survival'};

ReturnFn=@(l,aprime,a,rho,ej,b1,b2,b3,agej,J,r,delta,A,alpha,g_e)...
    KulishKentSmith2010_ReturnFn(l,aprime,a,rho,ej,b1,b2,b3,agej,J,r,delta,A,alpha,g_e);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

disp('Test ValueFnIter')
vfoptions=struct(); % Just use the defaults
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
toc

%% Distribution of agents
jequaloneDist=zeros([n_a,1],'gpuArray');
jequaloneDist(1,1)=1; % zero assets, not retired

% Distribution of agents in terms of ages
AgeWeightsParamNames={'mewj'};
Params.mewj=ones(N_j,1);
for jj=2:Params.J
    Params.mewj(jj)=Params.mewj(jj-1)/(1+Params.n);
end
for jj=Params.J+1:N_j
    Params.mewj(jj)=0;
end
Params.mewj=Params.mewj/sum(Params.mewj); % Normalize total mass to one
% KKS2010 never spell this out, but they seem to imply that this is what was used on page their 7.


%% Now solve the stationary dist to check things are set up correctly
disp('Test StationaryDist')
simoptions=struct(); % Just use the defaults
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
toc

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Functions to evaluate
FnsToEvaluate.K = @(l,aprime,a,agej,J) a*(agej<=J); % Aggregate assets (which is this periods state)
FnsToEvaluate.L = @(l,aprime,a,ej,agej,J) (1-l)*ej*(agej<=J); % Aggregate labour supply (in efficiency units)
% Note: because N_j is longer than actual lifespans, we need the (agej<=J) term

% General Equilibrium Equations
GEPriceParamNames={'r'};
GeneralEqmEqns.capitalmarket = @(r,L,K,A,alpha,delta) r-(A*alpha*(K^(alpha-1))*(L^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital

%% Test
disp('Test AggVars')
tic;
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);
toc

%% Solve for the initial general equilibrium
heteroagentoptions.verbose=1;
[p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
Params.r=p_eqm.r;

% Calculate some things about the equilibrium
[V_init, Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
StationaryDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

Params_init=Params;

save('./SavedOutput/KKS2010_InitialEqm.mat','V_init','Policy_init','StationaryDist_init','p_eqm','Params_init')

%% Set up parameters for the reform
T=150; % Number of periods to use when solving for transition path

for Reform=1:5 % See Table 3 of KKS2010, there are 4 reforms (they denote them 4.1 to 4.4, here they are just 1 to 4), 5th is needed for their Figure 10
    fprintf('Now doing Reform %i \n',Reform)
    Params=Params_init;
    % Note: not clear from paper what happens when J is increased to 70 (Reforms 3 and 4). I have
    % set it so there is initially noone ages 61-70, and then then
    % incrementally age into this. I guess this is what they did? But I
    % simply don't know what they did. (In principle this should be in
    % their code but I haven't gone looking).
    
    if Reform==1 % Baby Bust
        ParamPath.n=zeros(1,T); % We don't actually need the path on this
        ParamPath.mewj=zeros(N_j,T);
        % tt=1: shift all ages across one, and put in new first age
        ParamPath.mewj(1:Params.J,1)=[1+ParamPath.n(1);Params.mewj(1:Params.J-1)/Params.mewj(1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        for tt=2:T
            % Note that by using 1:Params.J it ensures there are always zero people aged over J
            ParamPath.mewj(1:Params.J,tt)=[1+ParamPath.n(tt);ParamPath.mewj(1:Params.J-1,tt-1)/ParamPath.mewj(1,tt-1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        end
        % Now normalize mewj for each period
        ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,1);
    elseif Reform==2 % Baby Boom and Bust
        % KKS2010, pg 15: "During the boom, agents act as if the higher fertility rate were to last forever. 
        %                  In other words, both of the changes in fertility are unanticipated."
        % This has to be modelled as two transition paths, put together.
        % Here we just set up the first, and later we do the second.
        ParamPath.n=0.024*ones(1,T);
        ParamPath.mewj=zeros(N_j,T);
        % tt=1
        ParamPath.mewj(1:Params.J,1)=[1+ParamPath.n(1);Params.mewj(1:Params.J-1)/Params.mewj(1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        for tt=2:T
            % Note that by using 1:Params.J it ensures there are always zero people aged over J
            ParamPath.mewj(1:Params.J,tt)=[1+ParamPath.n(tt);ParamPath.mewj(1:Params.J-1,tt-1)/ParamPath.mewj(1,tt-1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        end
        % Now normalize mewj for each period
        ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,1);
    elseif Reform==3 % Increase in longevity
        Params.J=70;
        Params.survival=[ones(1,Params.J),zeros(1,N_j-Params.J)]; % Note, this will just be all ones so is redundant but want to keep same setup as baseline
        % Need to set mewj to allow for those aged over 60
        ParamPath.mewj=Params.mewj.*ones(1,T);
        for tt=1:T
            for jjextra=1:min(tt,10) % 10 is the original J=60 to the new J=70
                ParamPath.mewj(60+jjextra,tt)=ParamPath.mewj(60+jjextra-1,tt)/(1+Params.n);
            end
        end
        % Now normalize mewj for each period
        ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,1);
        % Note: no need to recalculate v_jJ as this is already based on J (and b1,b2,b3), so is changing automatically
        % Note: no need to change 
    elseif Reform==4 % Reforms 2 and 3 together
        ParamPath.n=0.024*ones(1,T);
        Params.J=70; % Notice that I just use Params, rather than ParamPath (could do the later, but it is not necessary)
        Params.survival=[ones(1,Params.J),zeros(1,N_j-Params.J)]; % Note, this will just be all ones so is redundant but want to keep same setup as baseline
        ParamPath.mewj=zeros(N_j,T);
        % tt=1
        ParamPath.mewj(1:Params.J,1)=[1+ParamPath.n(1);Params.mewj(1:Params.J-1)/Params.mewj(1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        for tt=2:T
            % Note that by using 1:Params.J it ensures there are always zero people aged over J
            ParamPath.mewj(1:Params.J,tt)=[1+ParamPath.n(tt);ParamPath.mewj(1:Params.J-1,tt-1)/ParamPath.mewj(1,tt-1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        end
        % Now normalize mewj for each period
        ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,1);
        % Note: no need to recalculate v_jJ as this is already based on J (and b1,b2,b3)
    elseif Reform==5
        % This is not an actual reform. Is needed for Figure 10 which considers
        % what happens if survival is increased without improvement in 'health'
        Params.J=60; % Does not increase, this is controlling 'health' (via v_jJ)
        Params.survival=[ones(1,70),zeros(1,N_j-70)]; % This is the increase in survival to 70 periods
        % Need to set mewj to allow for those aged over 60
        ParamPath.mewj=Params.mewj.*ones(1,T);
        for tt=1:T
            for jjextra=1:min(tt,10) % 10 is the original 60 to the new 70
                ParamPath.mewj(60+jjextra,tt)=ParamPath.mewj(60+jjextra-1,tt)/(1+Params.n);
            end
        end
        % Now normalize mewj for each period
        ParamPath.mewj=ParamPath.mewj./sum(ParamPath.mewj,1);
    end
    
    
    %% Set up the general eqm conditions
    
    % The following transition path general equilibrium conditions for are just
    % the same as those the stationary model, but in some modeles they would differ
    GeneralEqmEqns_Transition.CapitalMarket = @(r,K,L,A,alpha,delta) r-(A*(alpha)*(K^(alpha-1))*(L^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
    
    transpathoptions.GEnewprice=3;
    % Need to explain to transpathoptions how to use the GeneralEqmEqns to
    % update the general eqm transition prices (in PricePath).
    transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
        {'CapitalMarket','r',0,0.1;...  % CaptialMarket GE condition will be positive if r is too big, so subtract
        };
    % Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
    % Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
    % A small 'factor' will make the convergence to solution take longer, but too large a value will make it
    % unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.
    
    transpathoptions.verbose=1;
    
    %% Solve the final general equilibrium
    Params_final=Params;
    if Reform==1
        Params_final.n=ParamPath.n(end);
        Params_final.mewj=ParamPath.mewj(:,end);
    elseif Reform==2 % Note: this is just the first of two stages
        Params_final.n=ParamPath.n(end);
        Params_final.mewj=ParamPath.mewj(:,end);
    elseif Reform==3
        Params_final.J=Params.J;
        Params_final.mewj=ParamPath.mewj(:,end);
        Params_final.survival=Params.survival;
    elseif Reform==4 % Note: this is just the first of two stages
        Params_final.n=ParamPath.n(end);
        Params_final.J=Params.J;
        Params_final.mewj=ParamPath.mewj(:,end);
        Params_final.survival=Params.survival;
    elseif Reform==5    
        Params_final.J=Params.J;
        Params_final.survival=Params.survival;
    end
    
    [p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_final, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    Params_final.r=p_eqm.r;
    
    % Calculate some things about the equilibrium
    [V_final, Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params_final, DiscountFactorParamNames, [],vfoptions);
    StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params_final,simoptions);
    
    %% Solve the transition path
    
    % Set up initial guess for the transition path
    PricePath0.r=[linspace(Params.r,Params_final.r,floor(T/2)),Params_final.r*ones(1,T-floor(T/2))];
    
    fprintf('Now computing the Transition path itself \n')
    PricePath=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, N_j, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightsParamNames, transpathoptions, simoptions, vfoptions);
    
    %% Now calculate some things about the transition path (The path for Value fn, Policy fn, and Agent Distribution)
    
    % You can calculate the value and policy functions for the transition path
    [VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions);
    
    % You can then use these to calculate the agent distribution for the transition path
    AgentDistPath=AgentDistOnTransPath_Case1_FHorz(StationaryDist_init,PricePath, ParamPath, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,pi_z,T, Params);
    
    %% If doing reform 2 or 4 this is actually just the first 20 years, and we then get a second surprise path.
    % KKS2010, pg 15: "During the boom, agents act as if the higher fertility rate were to last forever. 
    %                  In other words, both of the changes in fertility are unanticipated."
    
    if Reform==2
        Params.mewj=ParamPath.mewj(:,20); % Start from year 20 (need to use it to create ParamPath.mewj)
        StationaryDist_init=AgentDistPath(:,:,21); % The first period of the second stage of the reform
        
        % From year 21 onwards
        ParamPath2.n=zeros(1,T);
        ParamPath2.mewj=zeros(N_j,T);
        % tt=1
        ParamPath2.mewj(1:Params.J,1)=[1+ParamPath2.n(1);Params.mewj(1:Params.J-1)/Params.mewj(1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        for tt=2:T
            % Note that by using 1:Params.J it ensures there are always zero people aged over J
            ParamPath2.mewj(1:Params.J,tt)=[1+ParamPath2.n(tt);ParamPath2.mewj(1:Params.J-1,tt-1)/ParamPath2.mewj(1,tt-1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        end
        % Now normalize mewj for each period
        ParamPath2.mewj=ParamPath2.mewj./sum(ParamPath2.mewj,1);
        
        Params_final.n=ParamPath2.n(end);
        Params_final.mewj=ParamPath2.mewj(end);
        
    elseif Reform==4        
        Params.mewj=ParamPath.mewj(:,20); % Start from year 20 (need to use it to create ParamPath.mewj)
        StationaryDist_init=AgentDistPath(:,:,21); % The first period of the second stage of the reform
        
        % Note that J and survival just change once, it is only n (and hence mewj) that have the two stages.
        % J and survival are already in Params, so no need to repeat them here.
        
        % From year 21 onwards
        ParamPath2.n=zeros(1,T);
        ParamPath2.mewj=zeros(N_j,T);
        % tt=1
        ParamPath2.mewj(1:Params.J,1)=[1+ParamPath2.n(1);Params.mewj(1:Params.J-1)/Params.mewj(1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        for tt=2:T
            % Note that by using 1:Params.J it ensures there are always zero people aged over J
            ParamPath2.mewj(1:Params.J,tt)=[1+ParamPath2.n(tt);ParamPath2.mewj(1:Params.J-1,tt-1)/ParamPath2.mewj(1,tt-1)]; % Use that second age is normalized to one (as that way newborns are just 1*(1+n) )
        end
        % Now normalize mewj for each period
        ParamPath2.mewj=ParamPath2.mewj./sum(ParamPath2.mewj,1);
        
        Params_final.n=ParamPath2.n(end);
        Params_final.mewj=ParamPath2.mewj(end);
    end
    
    if Reform==2 || Reform==4
        % Solve the second stage of the transition, and then put the two transitions together to get what we are actually after.
                
        % Calculate the final stationary general equilibrium
        [p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_final, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
        Params_final.r=p_eqm.r;
        [V_final, Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params_final, DiscountFactorParamNames, [],vfoptions);
        StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params_final,simoptions);

        % Set up initial guess for the transition path
        PricePath0.r=[linspace(PricePath.r(21),Params_final.r,floor(T/2)),Params_final.r*ones(1,T-floor(T/2))];
        % Calculate the transition path (being lazy about updating PricePath0)
        PricePath2=TransitionPath_Case1_FHorz(PricePath0, ParamPath2, T, V_final, StationaryDist_init, n_d, n_a, n_z, N_j, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightsParamNames, transpathoptions, simoptions, vfoptions);
        
        % Need the agent distribution path as where it is in period 20 is our new initial agent distribution
        [VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions);
        AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(StationaryDist_init,PricePath, ParamPath, PolicyPath, AgeWeightsParamNames,n_d,n_a,n_z,N_j,pi_z,T, Params);

        % Put the two transitions together to get the actual path we are interested in
        ParamPath.n=[ParamPath.n(1:20), ParamPath2.n(1:T-20)];
        ParamPath.mewj=[ParamPath.mewj(:,1:20), ParamPath2.mewj(:,1:T-20)];
        size(ParamPath2.n)
        size(PricePath2.r)
        PricePath.r=[PricePath.r(1:20); PricePath2.r(1:T-20)];
        for tt=21:T
            AgentDistPath(:,:,tt)=AgentDistPath2(:,:,tt-20);
            VPath(:,:,tt)=VPath2(:,:,tt-20);
            PolicyPath(:,:,:,tt)=PolicyPath2(:,:,:,tt-20);
        end
    end
    
    
    %% Life-cycle profiles in initial and final equilibrium
    % Add some more FnsToEvaluate
    FnsToEvaluate.leisure = @(l,aprime,a,ej,agej,J) l*(agej<=J); % Leisure
    FnsToEvaluate.C = @(l,aprime,a,agej,J,ej,r,delta,A,alpha,g_e) (agej<=J)*KulishKentSmith2010_ConsumptionFn(l,aprime,a,ej,r,delta,A,alpha,g_e); % Consumption
    
    AgeConditionalStats_init=LifeCycleProfiles_FHorz_Case1(StationaryDist_init,Policy_init,FnsToEvaluate,[],Params_init,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeConditionalStats_final=LifeCycleProfiles_FHorz_Case1(StationaryDist_final,Policy_final,FnsToEvaluate,[],Params_final,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    
    save(['./SavedOutput/KKS2010_AgeCondStatsfinal_Reform',num2str(Reform),'.mat'],'AgeConditionalStats_init','AgeConditionalStats_final','Params_final','Params_init')
    
    if CreateFigures==1
        fig3=figure(3);
        subplot(3,1,1); plot(1:1:Params_init.J,AgeConditionalStats_init.K.Mean(1:Params_init.J),1:1:Params_final.J,AgeConditionalStats_final.K.Mean(1:Params_final.J))
        title('Life-Cycle Profile: Individual Wealth')
        legend('Initial','Final')
        subplot(3,1,2); plot(1:1:Params_init.J,AgeConditionalStats_init.leisure.Mean(1:Params_init.J),1:1:Params_final.J,AgeConditionalStats_final.leisure.Mean(1:Params_final.J))
        title('Life-Cycle Profile: Individual Leisure')
        subplot(3,1,3); plot(1:1:Params_init.J,AgeConditionalStats_init.C.Mean(1:Params_init.J),1:1:Params_final.J,AgeConditionalStats_final.C.Mean(1:Params_final.J))
        title('Life-Cycle Profile: Individual Consumption')

        if Reform==1
            saveas(fig3,'./SavedOutput/Graphs/KulishKentSmith2010_Figure5.png')
        elseif Reform==3
            saveas(fig3,'./SavedOutput/Graphs/KulishKentSmith2010_Figure8.png')
        end
    end
    
    
    %% Transition paths for some variables
    
    save(['./SavedOutput/KKS2010_Debug_Reform',num2str(Reform),'.mat'],'FnsToEvaluate','AgentDistPath','PolicyPath','PricePath','ParamPath','Params','Params_init','PricePath','Params', 'T', 'n_d', 'n_a', 'n_z', 'N_j', 'pi_z', 'd_grid', 'a_grid','z_grid', 'DiscountFactorParamNames', 'transpathoptions')    
    
    AggVars_init=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params_init, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);
    AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDistPath,PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, transpathoptions);
    % Just as a double check on the tails of the path
    AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params_final, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);
    
    
    KdivL_path=[AggVars_init.K.Mean; AggVarsPath.K.Mean]./[AggVars_init.L.Mean; AggVarsPath.L.Mean];
    w_path=Params.A*(1-Params.alpha)*(KdivL_path.^Params.alpha); % A and alpha dont change on path
    
    save(['./SavedOutput/KKS2010_KdivLratio_Reform',num2str(Reform),'.mat'],'KdivL_path','w_path','AggVars_init','AggVarsPath','Params','Params_init','PricePath','AggVars_final')

    if CreateFigures==1
        fig4=figure(4);
        subplot(3,1,1); plot(KdivL_path)
        legend('Capital to labor ratio')
        subplot(3,1,2); plot(w_path)
        legend('Real wage')
        subplot(3,1,3); plot([Params_init.r; PricePath.r])
        legend('Net real interest rate')
        title('Transitional Dynamics - Fall in Fertility')
        
        if Reform==1
            saveas(fig4,'./SavedOutput/Graphs/KulishKentSmith2010_Figure4.png')
        end
    end
    
    %% Growth rates of aggregates
    % We have the (per-capita) K and L, now inflate them to get the aggregates, then calculate their growth rates.
    
    if Reform==3
        ParamPath.n=Params.n*ones(1,T); % We didn't need a ParamPath.n for Reform 3, create one here just for the graph
    end
    % n is the growth rate of the model-age-one cohort
    % Need to calculate the growth rate of the agent mass from this
    nvec=[ParamPath.n(1); Params_init.n*ones(N_j-1,1)]; % Growth rates of each age-cohort
    PopulationGrowthRate=nan(1,T);
    
    for tt=1:T-1
        PopulationGrowthRate(tt)=sum(nvec.*ParamPath.mewj(:,tt));
        nvec=[ParamPath.n(tt+1); nvec(1:end-1)]; % Shift growth rates of age-cohorts by one, and put in growth rate of age-1 cohort
    end
    PopulationGrowthRate(T)=sum(nvec.*ParamPath.mewj(:,T));
    
    AggKpath=AggVarsPath.K.Mean.*cumprod(1+PopulationGrowthRate)';
    AggLpath=AggVarsPath.L.Mean.*cumprod(1+PopulationGrowthRate)';
    
    g_K=(AggKpath(2:end)-AggKpath(1:end-1))./AggKpath(1:end-1);
    g_L=(AggLpath(2:end)-AggLpath(1:end-1))./AggLpath(1:end-1);
    
    save(['./SavedOutput/KKS2010_GrowthRateAgg_Reform',num2str(Reform),'.mat'],'Params_init','g_K','g_L','T')

    if CreateFigures==1
        fig5=figure(5);
        plot(1:1:T,100*[Params_init.n;g_K],1:1:T,100*[Params_init.n;g_L]) % Note: in intitial stationary eqm all aggregates grow at rate n by definition
        legend('Growth rate of K_t','Growth rate of L_t')
        title('Fall in Fertility')

        if Reform==1
            saveas(fig5,'./SavedOutput/Graphs/KulishKentSmith2010_Figure6.png')
        end
    end
    
    %% Life-cycle profiles for those near retirement at time of change
    % Plot the life-cycle profiles for those aged 50, 55 and 60 at time of change.
    % Note that this will combine their pre-reform, with the transition path.
    
    % Calculate 'life-cycle profiles' over the transition
    AgeConditionalStatsPath=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPath,PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, simoptions);
    
    save(['./SavedOutput/KKS2010_LCP_Reform',num2str(Reform),'.mat'],'AgeConditionalStats_init','Params','AgeConditionalStatsPath')
    if Reform==3
        save(['./SavedOutput/KKS2010_LCP_Reform_extras',num2str(Reform),'.mat'],'AgeConditionalStats_init','Params','AgeConditionalStatsPath','FnsToEvaluate', 'AgentDistPath','PolicyPath', 'PricePath', 'ParamPath', 'Params', 'T', 'n_d', 'n_a', 'n_z', 'N_j', 'd_grid', 'a_grid','z_grid', 'simoptions')
    end
    
    if CreateFigures==1
        fig6=figure(6);
        % 'Leisure and Consumption Profiles with Increased Longevity for Agents Aged 50, 55 and 60 at Time of Change'
        subplot(3,1,1); plot(1:1:Params.J,[AgeConditionalStats_init.K.Mean(1:49),AgeConditionalStatsPath.K.initialages.Mean(50,50:Params.J)],'-b')
        hold on
        subplot(3,1,1); plot(1:1:Params.J,[AgeConditionalStats_init.K.Mean(1:54),AgeConditionalStatsPath.K.initialages.Mean(55,55:Params.J)],'-r')
        subplot(3,1,1); plot(1:1:Params.J,[AgeConditionalStats_init.K.Mean(1:59),AgeConditionalStatsPath.K.initialages.Mean(60,60:Params.J)],'-g')
        hold off
        title('Life-Cycle Profile: Individual Wealth')
        subplot(3,1,2); plot(1:1:Params.J,[AgeConditionalStats_init.leisure.Mean(1:49),AgeConditionalStatsPath.leisure.initialages.Mean(50,50:Params.J)],'-b')
        hold on
        subplot(3,1,2); plot(1:1:Params.J,[AgeConditionalStats_init.leisure.Mean(1:54),AgeConditionalStatsPath.leisure.initialages.Mean(55,55:Params.J)],'-r')
        subplot(3,1,2); plot(1:1:Params.J,[AgeConditionalStats_init.leisure.Mean(1:59),AgeConditionalStatsPath.leisure.initialages.Mean(60,60:Params.J)],'-g')
        hold off
        title('Life-Cycle Profile: Individual Leisure')
%         subplot(3,1,3); plot(1:1:50,AgeConditionalStats_init.C.Mean(1:50),'-b',50:1:Params.J,AgeConditionalStatsPath.C.initialages.Mean(50,50:Params.J),'-b')
        subplot(3,1,3); plot(1:1:Params.J,[AgeConditionalStats_init.C.Mean(1:49),AgeConditionalStatsPath.C.initialages.Mean(50,50:Params.J)],'-b')
        hold on
        subplot(3,1,3); plot(1:1:Params.J,[AgeConditionalStats_init.C.Mean(1:54),AgeConditionalStatsPath.C.initialages.Mean(55,55:Params.J)],'-r')
        subplot(3,1,3); plot(1:1:Params.J,[AgeConditionalStats_init.C.Mean(1:59),AgeConditionalStatsPath.C.initialages.Mean(60,60:Params.J)],'-g')
        legend('Age 50 at Reform','Age 55 at Reform','Age 60 at Reform','Location','northwest')
        hold off
        title('Life-Cycle Profile: Individual Consumption')
        
        if Reform==3
            saveas(fig6,'./SavedOutput/Graphs/KulishKentSmith2010_Figure9.png')
        end
    end
    
    
    

    %% Create Tables 4 and 5 (looping over solving the GE problem, but without the transition, just for Reform 1)
    if Reform==1
        Params=Params_init;
        nvec=[0.024,0.012,0,-0.012,-0.024];
        Jvec=[60,65,70,75,80];
        Table_CapitalIntensity=zeros(4,5);
        Table_ShareWorking=zeros(5,5);
        for n_c=1:length(nvec)
            Params.n=nvec(n_c);
            
            for J_c=1:length(Jvec)
                Params.J=Jvec(J_c);
            
                % For these Tables I use N_j=J (unlike in main codes)
                N_j_specific=Jvec(J_c);
                
                % Update the age-dependent parameters
                Params.agej=1:1:N_j_specific;
                Params.survival=[ones(1,Params.J),zeros(1,N_j_specific-Params.J)];
                if N_j_specific<=N_j
                    Params.ej=Params_init.ej(1:N_j_specific);
                else
                    Params.ej=[Params_init.ej,Params_init.ej(end)*ones(1,N_j_specific-N_j)];
                end
                
                % Changing n means we have to update mewj
                Params.mewj=ones(N_j_specific,1);
                for jj=2:Params.J
                    Params.mewj(jj)=Params.mewj(jj-1)/(1+Params.n);
                end
                Params.mewj=Params.mewj/sum(Params.mewj); % Normalize total mass to one
                % KKS2010 never spell this out, but they seem to imply that this is what was used on page their 7.

                
                [p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j_specific, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
                Params.r=p_eqm.r;
                % Calculate some things about the equilibrium
                [V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j_specific, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
                StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j_specific,pi_z,Params,simoptions);
                AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j_specific, d_grid, a_grid, z_grid);
                AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j_specific,d_grid,a_grid,z_grid,simoptions);
                
                if n_c<5
                    Table_CapitalIntensity(n_c,J_c)=AggVars.K.Mean/AggVars.L.Mean;
                end
                Table_ShareWorking(n_c,J_c)=100*sum(AgeConditionalStats.leisure.Mean)/Params.J;
            end
        end
        
        save('./SavedOutput/KulishKentSmith2010_Table_CapitalIntensity.mat','Table_CapitalIntensity')
        save('./SavedOutput/KulishKentSmith2010_Table_ShareWorking.mat','Table_ShareWorking')
        
        FID = fopen('./SavedOutput/LatexInputs/KulishKentSmith2010_Table4.tex', 'w');
        fprintf(FID, '\\begin{tabular}{lccccc} \\hline \n');
        fprintf(FID, 'n (percent) & \\multicolumn{5}{c}{T (years)} \\\\ \n');
        fprintf(FID, ' & %i & %i & %i & %i & %i \\\\ \\hline \n',Jvec(1),Jvec(2),Jvec(3),Jvec(4),Jvec(5));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(1), Table_CapitalIntensity(1,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(2), Table_CapitalIntensity(2,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(3), Table_CapitalIntensity(3,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(4), Table_CapitalIntensity(4,:));
        fprintf(FID, '\\hline \n \\end{tabular} \n');
        fclose(FID);
        
        FID = fopen('./SavedOutput/LatexInputs/KulishKentSmith2010_Table5.tex', 'w');
        fprintf(FID, '\\begin{tabular}{lccccc} \\hline \n');
        fprintf(FID, 'n (percent) & \\multicolumn{5}{c}{T (years)} \\\\ \n');
        fprintf(FID, ' & %i & %i & %i & %i & %i \\\\ \\hline \n',Jvec(1),Jvec(2),Jvec(3),Jvec(4),Jvec(5));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(1), Table_ShareWorking(1,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(2), Table_ShareWorking(2,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(3), Table_ShareWorking(3,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(4), Table_ShareWorking(4,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(5), Table_ShareWorking(5,:));
        fprintf(FID, '\\hline \n \\end{tabular} \n');
        fclose(FID);
        
        
        % Set Params back to what it was
        Params=Params_init;
    end
    
    %% Tables 6 and 7 require changing preferences and adding labor-augmenting productivity growth (just for Reform 1)
    if Reform==1
        
        Params=Params_init;
        Params.g_e=0.02; % Labor-augementing productivity growth
        Params.rho=1; % We need to switch to preferences that are compatible with a balanced-growth path
        
        % This is not obvious how to implement as paper is scarce on detail on this
        % extension (they appear to give enough to figure it out, just don't give actual equations).
        %
        % KKS2010 state: e_{j,t}=(1+g_e)e_{j,t-1}, where j is age and t is time
        % period, so that "each new generation faces a human capital profile during
        % their lifetime that is (1+g_e) times larger than the human capital
        % profile of the previous generation".
        %
        % My understanding involves adding a (1+g_e) term multiplying aprime in the
        % household problem, but leaving it otherwise unchanged. (log-log
        % preferences, which I personally hate, mean we just ignore the growth as
        % we can pull it out of the utility as an additive constant which as no
        % effect on policy, it changes the actual utility, but that is anyway only
        % defined up to a linear transformation and we are not doing any welfare
        % evaluation so who cares).
        %
        % If I understand correctly none of the aggregates (except consumption)
        % change expression, and nor is there any change to general equilibrium
        % conditions.
        
        nvec=[0.024,0.012,0,-0.012];
        Jvec=[60,65,70,75,80];
        Table6_CapitalIntensity=zeros(4,4);
        Table7_ShareWorking=zeros(4,4);
        for n_c=1:length(nvec)
            Params.n=nvec(n_c);
            
            for J_c=1:length(Jvec)
                Params.J=Jvec(J_c);
                N_j_specific=Jvec(J_c);
                
                % Update the age-dependent parameters
                Params.agej=1:1:N_j_specific;
                Params.survival=[ones(1,Params.J),zeros(1,N_j_specific-Params.J)];
                if N_j_specific<=N_j
                    Params.ej=Params_init.ej(1:N_j_specific);
                else
                    Params.ej=[Params_init.ej,Params_init.ej(end)*ones(1,N_j_specific-N_j)];
                end 
               
                % Changing n means we have to update mewj
                Params.mewj=ones(N_j_specific,1);
                for jj=2:Params.J
                    Params.mewj(jj)=Params.mewj(jj-1)/(1+Params.n);
                end
                Params.mewj=Params.mewj/sum(Params.mewj); % Normalize total mass to one
                % KKS2010 never spell this out, but they seem to imply that this is what was used on page their 7.
                
                [p_eqm,p_eqm_index, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j_specific, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
                Params.r=p_eqm.r;
                % Calculate some things about the equilibrium
                [V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j_specific, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
                StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j_specific,pi_z,Params,simoptions);
                AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j_specific, d_grid, a_grid, z_grid);
                AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j_specific,d_grid,a_grid,z_grid,simoptions);
                
                Table6_CapitalIntensity(n_c,J_c)=AggVars.K.Mean/AggVars.L.Mean;
                Table7_ShareWorking(n_c,J_c)=100*sum(AgeConditionalStats.leisure.Mean)/Params.J;
            end
        end
        
        save('./SavedOutput/KKS2010_Table6_CapitalIntensity.mat','Table6_CapitalIntensity')
        save('./SavedOutput/KKS2010_Table7_ShareWorking.mat','Table7_ShareWorking')
        
        FID = fopen('./SavedOutput/LatexInputs/KulishKentSmith2010_Table6.tex', 'w');
        fprintf(FID, '\\begin{tabular}{lccccc} \\hline \n');
        fprintf(FID, 'n (percent) & \\multicolumn{5}{c}{T (years)} \\\\ \n');
        fprintf(FID, ' & %i & %i & %i & %i & %i \\\\ \\hline \n',Jvec(1),Jvec(2),Jvec(3),Jvec(4),Jvec(5));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(1), Table6_CapitalIntensity(1,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(2), Table6_CapitalIntensity(2,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(3), Table6_CapitalIntensity(3,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(4), Table6_CapitalIntensity(4,:));
        fprintf(FID, '\\hline \n \\end{tabular} \n');
        fclose(FID);
        
        FID = fopen('./SavedOutput/LatexInputs/KulishKentSmith2010_Table7.tex', 'w');
        fprintf(FID, '\\begin{tabular}{lccccc} \\hline \n');
        fprintf(FID, 'n (percent) & \\multicolumn{5}{c}{T (years)} \\\\ \n');
        fprintf(FID, ' & %i & %i & %i & %i & %i \\\\ \\hline \n',Jvec(1),Jvec(2),Jvec(3),Jvec(4),Jvec(5));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(1), Table7_ShareWorking(1,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(2), Table7_ShareWorking(2,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(3), Table7_ShareWorking(3,:));
        fprintf(FID, '%8.1f &  %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*nvec(4), Table7_ShareWorking(4,:));
        fprintf(FID, '\\hline \n \\end{tabular} \n');
        fclose(FID);
        
        % Set Params back to what it was
        Params=Params_init;
    end
    
    
    
end % End of for loop over Reform=1:4






%% Figure 7 is based on all four reforms
if CreateFigures==1
    Figure7data=struct();
    for Reform=1:4
        load(['./SavedOutput/KKS2010_KdivLratio_Reform',num2str(Reform),'.mat'],'KdivL_path','AggVars_init','AggVars_final')
        Figure7data(Reform).KdivL_path=KdivL_path;
        % Reform 3 in Figure 7 goes the opposite direction to that in KKS2010.
        % As a double check, how does the capital-output ratio look in the final
        % and initial eqm for Reform 3?
        fprintf(['To check on Figure 7: for Reform ',num2str(Reform),' the initial capital-labor ratio was %8.3f and the final was %8.3f \n'],AggVars_init.K.Mean/AggVars_init.L.Mean ,AggVars_final.K.Mean/AggVars_final.L.Mean)
    end
    fig7=figure(7);
    plot(0:1:T,Figure7data(1).KdivL_path)
    hold on
    plot(0:1:T,Figure7data(2).KdivL_path)
    plot(0:1:T,Figure7data(3).KdivL_path)
    plot(0:1:T,Figure7data(4).KdivL_path)
    hold off
    legend('Reform 1','Reform 2','Reform 3','Reform 4')
    title('Capital-to-Labor Ratio')
    saveas(fig7,'./SavedOutput/Graphs/KulishKentSmith2010_Figure7.png')
end




%% Figure 10 is based on two reforms
if CreateFigures==1
    fig8=figure(8);
    Reform=5; load(['./SavedOutput/KKS2010_KdivLratio_Reform',num2str(Reform),'.mat'],'KdivL_path')
    subplot(3,1,1); plot(KdivL_path)
    title('Capital to labor ratio')
    Reform=1; load(['./SavedOutput/KKS2010_KdivLratio_Reform',num2str(Reform),'.mat'],'KdivL_path')
    hold on
    subplot(3,1,1); plot(KdivL_path)
    hold off
    legend('No health improvement','Health improvement')
    
    Reform=5; load(['./SavedOutput/KKS2010_AgeCondStatsfinal_Reform',num2str(Reform),'.mat'],'AgeConditionalStats_final','Params_final')
    subplot(3,1,2); plot(1:1:Params_final.J,AgeConditionalStats_final.leisure.Mean(1:Params_final.J))
    title('Individual leisure profile')
    Reform=1; load(['./SavedOutput/KKS2010_AgeCondStatsfinal_Reform',num2str(Reform),'.mat'],'AgeConditionalStats_final','Params_final')
    hold on
    subplot(3,1,2); plot(1:1:Params_final.J,AgeConditionalStats_final.leisure.Mean(1:Params_final.J))
    hold off
    
    Reform=5; load(['./SavedOutput/KKS2010_AgeCondStatsfinal_Reform',num2str(Reform),'.mat'],'AgeConditionalStats_final','Params_final')
    subplot(3,1,3); plot(1:1:Params_final.J,AgeConditionalStats_final.K.Mean(1:Params_final.J))
    title('Individual wealth profile')
    Reform=1; load(['./SavedOutput/KKS2010_AgeCondStatsfinal_Reform',num2str(Reform),'.mat'],'AgeConditionalStats_final','Params_final')
    hold on
    subplot(3,1,3); plot(1:1:Params_final.J,AgeConditionalStats_final.K.Mean(1:Params_final.J))
    hold off
    
    saveas(fig8,'./SavedOutput/Graphs/KulishKentSmith2010_Figure10.png')
end












