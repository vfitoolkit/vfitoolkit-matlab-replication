% Conesa & Krueger (1999) - Social Security Reform with Heterogeneous Agents

% I run this code twice, the first with the following vfoptions.fastOLG equal to one, and then a second time with this equal to zero.
vfoptions.fastOLG=0;

% The second time round I use the results of the first as an initial guess for the price paths:
if vfoptions.fastOLG==0
%     load ./SavedOutput/ConesaKrueger_FullResults_fastOLG.mat FullResults
    % I already computed the FullResults once through, so now just use that
    load ./SavedOutput/ConesaKrueger_FullResults.mat FullResults
end

% IndivProdShock=1; % 0 for 'no heterogeneity', 1 for Diaz-Gimenez et al (1997) the 'asymmetric case', 2 for Storesletten et al (1998) the 'symmetric case'.
%                   % Values of 3,4,5 are only used for Table 6 (2 state, 3 state, and 5 state) [In principle, 2 and 3 should be identical, but they differ slightly.]
%                   % Values of 6,7,8 are only used for Table 7 (2 states; 1 type, 2 type, 3 types, respectively)
%                   % Note that 6 is just the same as 3, but uses PType. Provides a double-check that PType is working correctly.
% PolicyReform=1; % 1=Policy Reform A,  2=Policy Reform B,  3=Policy Reform C

% Note: CK1999 model can actually be solved without decision variable. (You
% can calculate analytic value for labor supply based on a and aprime, see
% CK1999 paper. I model it explicitly so that code is easier to modify for
% other uses.)
if vfoptions.fastOLG==1
    n_d=51; % fraction of time worked
    n_a=301; % assets
else
    n_d=101; % fraction of time worked
    n_a=1001; % assets
end
% n_z % labour efficiency units, depends on 'IndivProdShock', so set below.
N_j=66; % age (number of periods, to be precise this is age-19)

% Note that equation (2.1) of CK1999 includes the conditional probability of survival in the expectations operator.
% In practice/codes it gets included in the discount factor. (This is clearer if you look at their equation (2.4)).

% CK1999 say they set tau=0.107 and b=0.5. Rather than follow this value for tau I treat
% it as a general equilibrium outcome (get a similar value in any case). As an initial 
% guess for tau I use that b=0.5 implies that tau=0.1189. As a double-check
% on the general equilibrium I just include it as part of the calculation
% of general equilibrium, and this can then be compared to the initial guess.
% (is a requirement of general equilibrium that can easily be
% calculated without actually solving model, see CK1999 Section 4 first
% paragraph. They verbally describe how you can relate tau and b from
% equations (2.12) and (2.13), this is implemented in my codes below when I
% set Params.tau_initial). [I also do this as for anything more than a
% flat tax rate this would have to be calculated as part of general eqm.]

% CK1999, when doing the 'symmetric' shocks case say that their Markov is
% approximating an AR(1) plus and iid. This is internally inconsistent with
% the actual model as if it were true then the agent's exogenous state
% should contain each of the AR(1) and the iid seperately as two different
% exogenous states (not just their sum eta). Further calculations suggest 
% that in fact this is not what happened, and CK1999 simply
% set eta in their model equal to exp(z) from equation (4.1) on page 767;
% another way to say this is that eta in equation (4.1) is a completely
% different eta from the eta in the model. Essentially, eta in eqn (4.1) is
% a typo; and exp(z) in eqn (4.1) is what in fact corresponds to eta in the
% model. These codes use exp(z) in eqn (4.1) as the process that is
% approximated to get eta. [Note: This is anyway irrelevant except of the
% results in Table 6.] [What CK1999 do, in ignoring epsilon and using their
% z from eqn (4.1) is standard practice based on interpretation of the iid
% epsilon as measurement error, it is simply notationally very confusing
% that eta is not eta!.] [CK1999 do not report Tauchen q, the
% hyperparameter for the Tauchen method, they may well in fact have used
% the Tauchen-Hussey method? (Tauchen-Hussey is bad idea: http://www.vfitoolkit.com/comment-on-guerrieri-lorenzoni-2017/ ]

% I do not have their original parameter values for s_j (which CK1999 called psi_j); but they use same source as Imorohoroglu, Imrohoroglu & Joines (1995) and I got IIJ1995 numbers so confident they are correct. 
% Not 100% sure on my epsilon_j, again I use IIJ1995 who used same source. 
% My tau is slightly different (see just above).
% I use tauchen method to discretize the shocks which is what they describe
% in paper, but their references section mentions Tauchen-Hussey so it is
% possible they used that instead.

% A few lines I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

% Set simoptions.ncores
simoptions.ncores=feature('numcores'); % Number of CPU cores

vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)
heteroagentoptions.verbose=1;

%% Set model parameters

%Preference parameters
Params.sigma=2;
Params.beta=0.97;
Params.gamma=0.42;

%Technology parameters
Params.alpha=0.36;
Params.delta=0.06;
Params.theta=1;

%People become economically active (j=1) at age 20, retire at 65, and max age is 85
Params.J=N_j; %=85-19
Params.Jr=46; % =65-19

%Population growth of 1.1%
Params.n=0.011;
%Social Security replacement rate
Params.b_initial=0.5; Params.b=Params.b_initial;

% Probability of surviving, conditional on being alive (Conesa & Krueger 1999 just say that they get them from Faber, 1982)
% This is same source as Imrohoroglu, Imrohoroglu & Joines (1995) and while replicating that paper they sent me a copy of the survival probabilities.
% The commented out lines below contain my original version using more recent numbers from 'same' source.
% I was sent the original survival probabilities by email by Selahattin
% Imrohoroglu and the contents of that 'psurv.dat' file are used below.
% %As this document is not available online I instead take the values from the updated version of the document, namely from 
% %F. C. Bell and M. L. Miller (2005), Life Tables for the United States Social Security Area 1900-2100, Actuarial Study No. 120, Office of the Chief Actuary
% %http://www.socialsecurity.gov/oact/NOTES/s2000s.html
% %Table 8 — Period Probabilities of Death Within One Year (qx) at Selected Exact Ages, by Sex and Calendar Year (Cont.)
% %The raw data from there is
% %          Sex and Exact Age
% %    |  Male                                                Female
% %Year| [0 30 60 65 70 100]                                  [0 30 60 65 70 100]
% %2010| [0.00587,0.00116,0.01086,0.01753,0.02785,0.39134]    [0.00495,0.00060,0.00734,0.01201,0.01912,0.34031]
% %I just take the numbers for Males, and then set my actual values based on a linear interpolation of the data.
% dj_temp=interp1([0,30,60,65,70,100],[0.00587,0.00116,0.01086,0.01753,0.02785,0.39134],0:1:100,'linear');
% Params.sj=1-dj_temp(20:85);
% Params.sj(1)=1;
% Params.sj(end)=0;
% % I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2005).
% % I have aditionally imposed that the prob of death at age 20 be zero and that prob of death at age 85 is one.
% Following are the orginal survival probabilites I was emailed by Selahattin Imrohoroglu in file 'psurv.dat'
Params.sj=[1.0000000, 0.99851000, 0.99844000, 0.99838000, 0.99832000, 0.99826000, 0.99820000, 0.99816000, 0.99815000, 0.99819000,...
    0.99826000, 0.99834000,0.99840000, 0.99843000, 0.99841000, 0.99835000, 0.99828000, 0.99818000, 0.99807000, 0.99794000,...
    0.99778000, 0.99759000, 0.99737000, 0.99712000, 0.99684000, 0.99653000, 0.99619000, 0.99580000, 0.99535000, 0.99481000,...
    0.99419000, 0.99350000, 0.99278000, 0.99209000, 0.99148000, 0.99088000, 0.99021000, 0.98942000, 0.98851000, 0.98746000,...
    0.98625000, 0.98495000, 0.98350000, 0.98178000, 0.97974000, 0.97743000, 0.97489000, 0.97226000, 0.96965000, 0.96715000,...
    0.96466000, 0.96200000, 0.95907000, 0.95590000, 0.95246000, 0.94872000, 0.94460000, 0.94017000, 0.93555000, 0.93077000, ...
    0.92570000, 0.92030000, 0.91431000, 0.90742000, 0.89948000];
Params.sj=[Params.sj,0]; % CK1999 has one more period than IIJ1995, but this period anyway ends in certain death so I am not 'missing' a conditional surivival probability for the final age

% Age profile of productivity (based on lifetime profile of earnings, Hansen 1993), Epsilon_j
% I have to interpolate this in some manner to get the values for each age. I am not sure how CK1999 did the interpolation as they never mention it.
% Again this is same source as Imrohoroglu, Imrohoroglu & Joines (1995) and
% while replicating that paper they sent me a copy of the their interpolated numbers which I will assume are same as CK1999.
%First the raw data of Hansen 1993:
%Table II. Weights assigned to age-sex
%         Males       Females
%Age      Weight      Weight
%16-19    0.56        0.52
%20-24    0.78        0.69
%25-34    1.14        0.89
%35-44    1.37        0.90
%45-54    1.39        0.87
%55-64    1.33        0.84
%65 +     0.89        0.66
%14-17    0.56        0.52
%14-19    0.56        0.52
%18-24    0.78        0.69
%25-44    1.24        0.89
%45-64    1.37        0.86
%epsilon_j is set based on the first entries of the data for males
% The following commented out part is my original attempt, below that are a
% copy-paste of the numbers used by IIJ1995 which I was sent by email by
% Selahattin Imrohoroglu.
% epsilon_j=zeros(Params.J,1); 
% epsilon_j(1:5)=0.78; epsilon_j(6:15)=1.14; epsilon_j(16:25)=1.37; 
% epsilon_j(26:35)=1.39; epsilon_j(36:45)=1.33; epsilon_j(46:66)=0.89;
% Params.epsilon_j=epsilon_j/epsilon_j(1); %Conesa & Krueger refer to labour efficiency as 1 at age 20, so presumably they make some renormalization like this
% The following line is the contents of 'comboeff.dat' file I received from Selahattin Imrohoroglu (see a few lines above)
Params.epsilon_j=[0.36928407;  0.42120465;  0.49398951; 0.56677436; 0.63955922; 0.71234407; 0.78512892; 0.81551713; 0.84590532; 0.87629354; 0.90668176; 0.93706995; 0.96379826; 0.99052656; 1.0172549; 1.0439831; 1.0707114; 1.0819067; 1.0931018; 1.1042970; 1.1154923; 1.1266874; 1.1435058; 1.1603241; 1.1771424; 1.1939607; 1.2107789; 1.2015913; 1.1924037; 1.1832160; 1.1740283; 1.1648407; 1.1609988; 1.1571568; 1.1533149; 1.1494730; 1.1456310; 1.1202085; 1.0947859; 1.0693634; 1.0439408; 1.0185183; 0.99309576; 0.96767321];
% CK1999 have one more time/age period than IIJ1995 did, so I will repeat the last entry twice
Params.epsilon_j=[Params.epsilon_j;Params.epsilon_j(end)];
Params.epsilon_j=[Params.epsilon_j; zeros(Params.J-Params.Jr+1,1)]; % Fill in the zeros for retirement age
Params.epsilon_j=Params.epsilon_j/Params.epsilon_j(1); %Conesa & Krueger refer to labour efficiency as 1 at age 20, so presumably they make some renormalization like this
% Conesa & Krueger ban retirees from working. One way to do this is just to
% set epsilon_j=0 for retirees. Rather than do this via epsilon_j it is
% done by I_j which in practice is multiplied by epsilon_j.
Params.I_j=[ones(Params.Jr-1,1);zeros(Params.J-Params.Jr+1,1)];
% Note that with my original version I_j was needed, but is redundant using the original epsilon_j numbers from IIJ1995


vfoptions.verbose=0;
vfoptions.policy_forceintegertype=2; % Policy was not being treated as integers (one of the elements was 10^(-15) different from an integer)
vfoptions.lowmemory=0; % This is changed to other values for some of the model variants.
% simoptions.nsims=4*10^5;
% simoptions.iterate=1;
heteroagentoptions.verbose=1;

%%
for IndivProdShock=6:8 % 0:8
    
    %Idiosycratic productivity shocks
    if IndivProdShock<=5
        Names_i={}; % Not used, this is needed when want different agent types as in IndivProdShock 6 to 8.
        PTypeDistParamNames={}; % Not used, this is needed when want different agent types as in IndivProdShock 6 to 8.
        
        % Parameters in Table IV: 1 for Storesletten et al (1998) the 'symmetric case', 2 for Diaz-Gimenez et al (1997) the 'asymmetric case'.
        if IndivProdShock==0 % 0 for the 'no heterogeneity'
            % CK1999 never explicitly describe what 'no hetergeneity' means. But
            % given that idiosyncratic productivity shocks eta are the source of
            % (ex-post) heterogeneity presumably it is simply a matter of setting
            % them to their mean value of 1.
            eta_grid=1;
            pi_eta=1;
        elseif IndivProdShock==1 % 1 for Storesletten et al (1998) the 'symmetric case'.
            eta_grid=[0.73; 1.27];
            pi_eta=[0.82, 0.18; 0.18, 0.82];
            % I should get the same using Tauchen method
            % Following lines are needed just to create Table 4
            Params.eta1=eta_grid(1); Params.eta2=eta_grid(2);
            Params.pi1=pi_eta(1,1); Params.pi2=pi_eta(2,2);
        elseif IndivProdShock==2 % 2 for Diaz-Gimenez et al (1997) the 'asymmetric case'
            eta_grid=[0.5;3];
            pi_eta=[0.9811, 1-0.9811; 1-0.9261, 0.9261];
            % Following lines are needed just to create Table 4
            Params.eta1=eta_grid(1); Params.eta2=eta_grid(2);
            Params.pi1=pi_eta(1,1); Params.pi2=pi_eta(2,2);
        elseif IndivProdShock>=3 && IndivProdShock<=5 % 3 'symmetric case', Tauchen approx of AR(1) from Storesletten et al (1998)
            Params.Tauchen_q=0.45; % For the two state this seems to do the trick. I will assume this is kept constant as the number of grid points in changed.
            % Note to self: I suspect that while paper says they use
            % Tauchen method they may in fact have used the Tauchen-Hussey
            % method as it is the later paper that appears in the
            % refererences?
            tauchenoptions.parallel=1; % use cpu rather than gpu for Tauchen method
            % The following three parameters are just used to determine eta_grid and pi_eta, both (eta_grid and pi_eta are themselves given directly in paper, pg 767-8 for the main 'symmetric case' hence they were not needed there).
            Params.rho=0.935;
            Params.sigmasq_epsilon=0.0017;
            Params.sigmasq_upsilon=0.061;
            if IndivProdShock==3
                % [This should be the same as the IndivProdShock==2 case.]
                n_z=2;
            elseif IndivProdShock==4
                n_z=3;
            elseif IndivProdShock==5
                if vfoptions.fastOLG==1
                    vfoptions.lowmemory=1;
                end
                n_z=5;
            end
            % What CK1999 call z I have called u1
            [u_grid, pi_u]=TauchenMethod(0,Params.sigmasq_upsilon,Params.rho,n_z,Params.Tauchen_q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,[tauchenoptions]), transmatix is (z,zprime)
            %     [epsilon_grid, pi_epsilon]=TauchenMethod(0,Params.sigmasq_epsilon,0,n_z,q); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,[tauchenoptions]), transmatix is (z,zprime)
            %     u_grid=u1_grid+epsilon_grid;
            %     eta_grid=exp(u_grid);
            eta_grid=exp(u_grid);
            pi_eta=pi_u;
        end
        n_z=length(eta_grid);
        z_grid=eta_grid;
        pi_z=pi_eta;
    end
    
    % There are also some where fixed types are used for Table 7.
    % The VFI Toolkit has a 'Permanent Type' designed for allowing for
    % genuinely different agent types.
    if IndivProdShock>=6
        % See CK1999, page 789 for explanation.
        % Same issue around notation and eta in model actually corresponding to the (exponential of) z in equation (6.1) [Not the eta in eqn (6.1)]
        Params.Tauchen_q=0.45; % For the two state this seems to do the trick.
        tauchenoptions.parallel=1; % use cpu rather than gpu for Tauchen method
        % The following three parameters are just used to determine eta_grid and pi_eta, both (eta_grid and pi_eta are themselves given directly in paper, pg 767-8 for the main 'symmetric case' hence they were not needed there).
        Params.rho=0.98;
        Params.sigmasq_epsilon=0.005; % This is never actually used for anything
        Params.sigmasq_upsilon=0.019;
        Params.sigmasq_alpha=0.326;
        n_z=2;
        % What CK1999 call z I have called u
        [u_grid, pi_u]=TauchenMethod(0,Params.sigmasq_upsilon,Params.rho,n_z,Params.Tauchen_q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,[tauchenoptions]), transmatix is (z,zprime)
%         eta_grid=exp(u_grid);
        pi_eta=pi_u;
        z_grid=u_grid; % Note that for Ptype z_grid uses u_grid rather than eta_grid (eta=exp(u))
        pi_z=pi_eta; % Note: these appear in 'reverse' order.
        % CK1999 call the permanent types 'alpha', they also use 'alpha' as
        % a parameter of the production fn. For this reason I call the permanent types 'zeta'.
        PTypeDistParamNames={'PTypeWeights'};
        if IndivProdShock==6 % One fixed type: note that this is actually almost exactly the the same things as IndivProdShock==3 (minor difference in Params.rho and Params.sigmadq_upsilon); I doesn't need to be done using PType but I do to make it eaiser to compare to IndivProdShock==7 and 8
            Names_i={'type1'};
            Params.zeta.type1=0;
            Params.PTypeWeights=[1]; % All agents are of the one type
        elseif IndivProdShock==7 % Two fixed types
            Names_i={'type1','type2'};
            Params.zeta.type1=-0.43*sqrt(Params.sigmasq_alpha);
            Params.zeta.type2=0.43*sqrt(Params.sigmasq_alpha);
            Params.PTypeWeights=[0.5; 0.5]; % Make agents half of each type
        elseif IndivProdShock==8 % Three fixed types
            Names_i={'type1','type2','type3'};
            Params.zeta.type1=-0.675*sqrt(Params.sigmasq_alpha);
            Params.zeta.type2=0*sqrt(Params.sigmasq_alpha);
            Params.zeta.type3=0.675*sqrt(Params.sigmasq_alpha);
            Params.PTypeWeights=[1/3; 1/3; 1/3]; % Make agents one-third of each type
        end
        
    end
    N_z=prod(n_z);

    %% Grid for assets
    % k_grid=55*(linspace(0,102,n_a).^3)';
    % When n_a=102 this is the actual grid chosen by Conesa & Krueger (1999) according to their description (pg 793). Note though that they use linear interpolation for the decision variable.
    % I suspect that their description is erroneous and is meant to be
    k_grid=55*(linspace(0,1,n_a).^3)';
    % Otherwise, strictly following their description, the maximum value of k_grid is so huge compared to maximum (period) income it seems silly.
    
    %% Get problem into format for using the toolkit
    d_grid=linspace(0,1,n_d)'; %fraction of time worked
    a_grid=k_grid; %assets

    %% Calculate the population distribution across the ages (population at time
    % t is made stationary by dividing it by \Pi_{i=1}^{t} (1+n_{i}) (product of all
    % growth rates up till current period))
    % CK2009 do not appear to give it an name (they instead say in
    % Computational Appendix that they simply work with the model after making
    % all the adjustments for population growth). The VFI Toolkit needs to give
    % it a name so that it can be automatically used when calculating model
    % outputs.
    Params.mewj=ones(Params.J,1);
    for jj=2:Params.J
        Params.mewj(jj)=Params.mewj(jj-1)*(1/(1+Params.n))*Params.sj(jj-1);
    end
    Params.mewj=Params.mewj/sum(Params.mewj); % normalize to measure one
    Params.mewj=Params.mewj'; % Age weights must be a row vector.
    
    AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.
    
    % This population distribution across the ages is also used to back out initial guess for the calibrated value for tau.
    Params.tau_initial=Params.b*sum(Params.mewj(:,Params.Jr:Params.J))/sum(Params.mewj(:,1:Params.Jr-1)); % From eqn 2.13 and 2.12 we get: tau=b*(frac of population retired)/(fraction of population of working age)
    % Note that rather than use this directly as the value for tau I instead
    % compute it as a general equilbrium condition. This is done just as a double check on rounding errors.
    
    % Stationary distribution of eta (the exogenous state z, the efficiency labour process)
    statdist_z=((ones(1,N_z)/N_z)*(pi_z^(10^6)))'; % Not really needed for anything. Just calculate it out of interest.

    %% General eqm variables: give some initial values
    GEPriceParamNames={'r','tau','SS','Tr'};
    Params.r=0.06; % interest rate on assets
    Params.tau=Params.tau_initial; % Payroll tax rate.
    Params.SS=0.4; % Benefits level for retirees
    Params.Tr=0.3; % lumpsum transfers made out of the accidental bequests
    % These were originally r=0.06, SS=1.2, Tr=1. Changed to present values as
    % these are much closer to the solution and so substantially reduce the run
    % time for finding General Eqm.

    %% Now, create the return function
    DiscountFactorParamNames={'beta','sj'};
    
    if IndivProdShock<=5
        ReturnFn=@(l,aprime,a,eta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma) ConesaKrueger1999_ReturnFn(l,aprime,a,eta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma)
        ReturnFnParamNames={'r','tau','epsilon_j','I_j','Tr','SS','theta','alpha','delta','gamma','sigma'}; %It is important that these are in same order as they appear in 'ConesaKrueger1999_ReturnFn'
    elseif IndivProdShock>=6
        % Uses u & zeta as the exogenous shocks instead of eta (note that eta=exp(u))
        ReturnFn=@(l,aprime,a,u,zeta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma) ConesaKrueger1999_Ptype_ReturnFn(l,aprime,a,u,zeta,r,tau,epsilon_j,I_j,Tr,SS,theta,alpha,delta,gamma,sigma)
        ReturnFnParamNames={'zeta','r','tau','epsilon_j','I_j','Tr','SS','theta','alpha','delta','gamma','sigma'}; %It is important that these are in same order as they appear in 'ConesaKrueger1999_ReturnFn'
    end
    %% Initial distribution of agents at birth (j=1)
    jequaloneDist=zeros([n_a,n_z]);
    jequaloneDist(1,:)=statdist_z; % All agents born with zero assets and with based on stationary distribution of the exogenous process on labour productivity units (comments elsewhere in CK1999 suggest that this is what they did, is not explicit in the paper)

    %% Some things for final general equilibrium
    % All three policy reforms have the same end point, namely completely terminate the social security systen:
    Params.b_final=0;
    % We know what the General Eqm values for tau and SS will be, so may as well set this as our initial guess.
    Params.tau_final=0;
    Params.SS_final=0;
    
    %% The transition paths for the policy reforms

    % Number of time periods to allow for the transition (if you set T too low
    % it will cause problems, too high just means run-time will be longer).
    % (Problems in the sense that solution would be 'incorrect', it would likely still solve)
    T=150; % CK1999 use T=150
    
    % Because we are going to calculate the general equilibrium transition
    % paths lets set it so that the general eqm (initial and final) are
    % calculated to very high (likely excessive) levels of accuracy.
    heteroagentoptions.toleranceGEprices=10^(-5); % Final eqm prices need to be highly accurate when using transition paths
    heteroagentoptions.toleranceGEcondns=10^(-5); % Final eqm condns need to be highly accurate when using transition paths
    
    for PolicyReform=1:3
        % To indicate progress:
        save ./SavedOutput/ConesaKrueger1999_ProgressIndicator.mat IndivProdShock PolicyReform
        
        % Three reforms are considered. These represent three different paths for the social security benefits (replacement rate) b.
        % ParamPathNames={'b'}; % This is the parameter that gets changed 'away' from it's initial value.
        % Note: The following quotes from CK1999 use a different timing convention for periods. Their period 2 is period 1 in the VFI Toolkit timing.
        if PolicyReform==1 % 1=Policy Reform A
            % "Beginning with period 2 the replacement rate is set equal to 0 and
            % stays there forever, i.e., b_1=0.5, b_t=0 for all t>1. This reform
            % terminates the social security system immediately and does not honor
            % entitlements to social security payments."
            ParamPath.b=zeros(T,1); % ParamPath is matrix of size T-by-'number of parameters that change over path'
        elseif PolicyReform==2 % 2=Policy Reform B
            % "The replacement rate, 50%, is linearly reduced by one percentage
            % point over 50 periods, i.e., bt=0.5-0.01(t-1), t=1,2,...,50, bt=0 for
            % all t>50. This reform terminates the social security system gradually
            % so that entitlements are partially honored and payroll taxes are
            % accordingly reduced to finance progressively smaller benefits."
            ParamPath.b=[linspace(0.49,0,50)'; zeros(T-50,1)]; % ParamPath is matrix of size T-by-'number of parameters that change over path'
            % Note: Start at 0.49 as the "50 periods" of CK1999 includes the initial steady state prior to announcement.
        elseif PolicyReform==3 %  3=Policy Reform C
            % "The replacement rate is fixed for 20 years at 50% and set to 0
            % thereafter, i.e., bt=0.5, t=1,2,...,20, bt=0 for all t>20. Therefore,
            % all individuals retired or about to retire keep their social security
            % benefits, but future retirees anticipate that they will receive only
            % part or no social security benefits. This reform allows agents to
            % readjust their plans for the anticipated reform in 20 years."
            ParamPath.b=[0.5*ones(19,1); zeros(T-19,1)]; % ParamPath is matrix of size T-by-'number of parameters that change over path'
            % Note: 19 as the "20 periods" of CK1999 includes the initial steady state prior to announcement.
        end
        
        % We need to give an initial guess for the price path on interest rates.
        if vfoptions.fastOLG==0
            PricePath0=FullResults(IndivProdShock+1,PolicyReform).Output.PricePath; % This is the version that was calculated with vfoptions.fastOLG==1
        else
            PricePath0=struct(); % Is filled in internally by ConesaKrueger1999_Fn.m based on the initial and final eqm
        end
        % PricePathNames={'r','tau','SS','Tr'};
        
        % Now just run the TransitionPath_Case1 command (all of the other inputs
        % are things we had already had to define to be able to solve for the
        % initial and final equilibria)
        transpathoptions.weightscheme=3; % default =1
        transpathoptions.verbose=1;
%         vfoptions.fastOLG=1; % Set at top of this script
        transpathoptions.maxiterations=200; % default is ???
        transpathoptions.oldpathweight=0.8; % default =0.9
        transpathoptions.historyofpricepath=1;
        vfoptions.policy_forceintegertype=2;


        transpathoptions.GEnewprice=1;
        % Setting transpathoptions.GEnewprice to 1 means the GeneralEqmEqns have instead been expressed as formulae to create 'new updated general equilibrium variables'.
        % This is much faster for the convergence of the transition path, which is otherwise quite complex.
        % Notice that in the present example this is very easy as each of the
        % original GeneralEqmEqns was just one of theses prices minus 'something',
        % so we now just rewrite them as 'something'.
        GeneralEqmEqnParamNames=struct();
        GeneralEqmEqnParamNames(1).Names={'theta','alpha','delta'};
        GeneralEqmEqn_1 = @(AggVars,GEprices,theta,alpha,delta) GEprices(1)-0.5*(GEprices(1)-(theta*(alpha)*(AggVars(1)^(alpha-1))*(AggVars(2)^(1-alpha))-delta)); % Rate of return on assets is related to Marginal Product of Capital
        GeneralEqmEqnParamNames(2).Names={'b'};
        GeneralEqmEqn_2 = @(AggVars,GEprices,b) GEprices(2)-1*(GEprices(2)-b*(1-AggVars(4))/AggVars(4)); % From eqn 2.13 and 2.12 we get: tau=b*(frac of population retired)/(fraction of population of working age)
        GeneralEqmEqnParamNames(3).Names={'b','theta','alpha'};
        GeneralEqmEqn_3 = @(AggVars,GEprices,b,theta,alpha) GEprices(3)-0.5*(GEprices(3)-b*(theta*(1-alpha)*((AggVars(1)/AggVars(2))^alpha))*AggVars(2)/AggVars(4)); % Social Security adds up, based on eqn (2.12) b*w*N/(fraction working age) [Note that because of how tau is calibrated (from b) this also means eqn (2.13) will hold.]
        GeneralEqmEqnParamNames(4).Names={};
        GeneralEqmEqn_4 = @(AggVars,GEprices) GEprices(4)-0.5*(GEprices(4)-AggVars(3)); % Accidental bequests (adjusted for population growth) are equal to transfers received
        GeneralEqmEqns_Transition={GeneralEqmEqn_1, GeneralEqmEqn_2, GeneralEqmEqn_3, GeneralEqmEqn_4};
        
        FullResults(IndivProdShock+1,PolicyReform).Output=ConesaKrueger1999_Fn(IndivProdShock,Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid,pi_z,ReturnFn, ReturnFnParamNames, DiscountFactorParamNames, jequaloneDist, GEPriceParamNames, AgeWeightsParamNames, T, ParamPath, PricePath0, GeneralEqmEqns_Transition, vfoptions, simoptions, heteroagentoptions, transpathoptions, Names_i, PTypeDistParamNames);
        % Notice the IndivProdShock+1
        
        ReplicationStatus=gather([IndivProdShock, PolicyReform]);
        save ./SavedOutput/ConesaKrueger_ReplicationStatus.mat ReplicationStatus
        if vfoptions.fastOLG==0
            save ./SavedOutput/ConesaKrueger_FullResults.mat FullResults a_grid -v7.3
        else
            save ./SavedOutput/ConesaKrueger_FullResults_fastOLG.mat FullResults a_grid -v7.3
        end
    end
end

%%
% load ./SavedOutput/ConesaKrueger_FullResults_fastOLG.mat FullResults
load ./SavedOutput/ConesaKrueger_FullResults.mat FullResults a_grid
T=150

%% Now to create all the Tables and Figures themselves.

%% Tables 1 to 4 are just the parameter values of the model. I still replicate these but it is somewhat pointless to do so.

% Table 1
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table1.tex', 'w');
fprintf(FID, 'Preference Parameters \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}cc} \n \\hline \n');
fprintf(FID, ' Parameter & Value  \\\\ \\hline \n');
fprintf(FID, ' $\\sigma$ & %8.2f \\\\ \n', Params.sigma);
fprintf(FID, ' $\\beta$ & %8.2f \\\\ \n', Params.beta);
fprintf(FID, ' $\\gamma$ & %8.2f \\\\ \n', Params.gamma);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fclose(FID);

% Table 2
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table2.tex', 'w');
fprintf(FID, 'Demographics \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}cc} \n \\hline \n');
fprintf(FID, ' Parameter & Value  \\\\ \\hline \n');
fprintf(FID, ' $J$ & %d \\\\ \n', Params.J);
fprintf(FID, ' $j_r$ & %d \\\\ \n', Params.Jr);
fprintf(FID, ' $s_j$ & Bell and Miller (2005) \\\\ \n');
fprintf(FID, ' $n$ & %8.3f \\\\ \n', Params.n);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Renamed the age-conditional surivial probability as $s_j$, while CK1999 call it $\\psi_j$. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 3
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table3.tex', 'w');
fprintf(FID, 'Technology Parameters \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}cc} \n \\hline \n');
fprintf(FID, ' Parameter & Value  \\\\ \\hline \n');
fprintf(FID, ' $\\alpha$ & %8.2f \\\\ \n', Params.alpha);
fprintf(FID, ' $\\delta$ & %8.2f \\\\ \n', Params.delta);
fprintf(FID, ' $\\theta$ & %d \\\\ \n', Params.theta);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fclose(FID);

% Table 4
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table4.tex', 'w');
fprintf(FID, 'Individual Productivity Parameters \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}ccc} \n \\hline \n');
fprintf(FID, ' Parameter & Diaz-Jimenez et al. & Storeslettern et al.  \\\\ \\hline \n');
fprintf(FID, ' $\\eta_1$ & %8.1f & %8.2f \\\\ \n', FullResults(3,1).Output.Params.eta1, FullResults(2,1).Output.Params.eta1);
fprintf(FID, ' $\\eta_2$ & %8.1f & %8.2f \\\\ \n', FullResults(3,1).Output.Params.eta2, FullResults(2,1).Output.Params.eta2);
fprintf(FID, ' $\\pi_1$ & %8.4f & %8.2f \\\\ \n', FullResults(3,1).Output.Params.pi1, FullResults(2,1).Output.Params.pi1);
fprintf(FID, ' $\\pi_2$ & %8.4f & %8.2f \\\\ \n', FullResults(3,1).Output.Params.pi2, FullResults(2,1).Output.Params.pi2);
fprintf(FID, ' $\\epsilon_j$ & Hansen (1993) & Hansen (1993) \\\\ \n');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fclose(FID);

%% Table 5

%Table 5
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table5.tex', 'w');
fprintf(FID, 'Steady-state Results \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{2}{c}{No heterogeneity} & \\multicolumn{2}{c}{Het. (sym. case)} & \\multicolumn{2}{c}{Het. (asym. case)} \\\\ \\cline{2-3} \\cline{4-5} \\cline{6-7} \n');
fprintf(FID, 'Variables & Init. St. St. & Fin. St. St. & Init. St. St. & Fin. St. St. & Init. St. St. & Fin. St. St. \\\\ \\hline \n');
fprintf(FID, 'b   & %i \\%% & %i \\%% & %i \\%% & %i \\%% & %i \\%% & %i \\%% \\\\ \n', 100*FullResults(1,1).Output.StationaryEqmStats(1).b, 100*FullResults(1,1).Output.StationaryEqmStats(2).b, 100*FullResults(2,1).Output.StationaryEqmStats(1).b, 100*FullResults(2,1).Output.StationaryEqmStats(2).b, 100*FullResults(3,1).Output.StationaryEqmStats(1).b, 100*FullResults(3,1).Output.StationaryEqmStats(2).b);
fprintf(FID, 'r   & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% \\\\ \n', 100*FullResults(1,1).Output.StationaryEqmStats(1).r, 100*FullResults(1,1).Output.StationaryEqmStats(2).r, 100*FullResults(2,1).Output.StationaryEqmStats(1).r, 100*FullResults(2,1).Output.StationaryEqmStats(2).r, 100*FullResults(3,1).Output.StationaryEqmStats(1).r, 100*FullResults(3,1).Output.StationaryEqmStats(2).r);
fprintf(FID, 'w   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\%% \\\\ \n', FullResults(1,1).Output.StationaryEqmStats(1).w, FullResults(1,1).Output.StationaryEqmStats(2).w, FullResults(2,1).Output.StationaryEqmStats(1).w, FullResults(2,1).Output.StationaryEqmStats(2).w, FullResults(3,1).Output.StationaryEqmStats(1).w, FullResults(3,1).Output.StationaryEqmStats(2).w);
fprintf(FID, 'h   & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% & %8.1f \\%% \\\\ \n', FullResults(1,1).Output.StationaryEqmStats(1).h, FullResults(1,1).Output.StationaryEqmStats(2).h, FullResults(2,1).Output.StationaryEqmStats(1).h, FullResults(2,1).Output.StationaryEqmStats(2).h, FullResults(3,1).Output.StationaryEqmStats(1).h, FullResults(3,1).Output.StationaryEqmStats(2).h);
fprintf(FID, 'K/Y & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(1,1).Output.StationaryEqmStats(1).Kdivy, FullResults(1,1).Output.StationaryEqmStats(2).Kdivy, FullResults(2,1).Output.StationaryEqmStats(1).Kdivy, FullResults(2,1).Output.StationaryEqmStats(2).Kdivy, FullResults(3,1).Output.StationaryEqmStats(1).Kdivy, FullResults(3,1).Output.StationaryEqmStats(2).Kdivy);
fprintf(FID, 'y   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(1,1).Output.StationaryEqmStats(1).y, FullResults(1,1).Output.StationaryEqmStats(2).y, FullResults(2,1).Output.StationaryEqmStats(1).y, FullResults(2,1).Output.StationaryEqmStats(2).y, FullResults(3,1).Output.StationaryEqmStats(1).y, FullResults(3,1).Output.StationaryEqmStats(2).y);
fprintf(FID, 'SS/y   & %8.1f \\%% & %8.0f \\%% & %8.1f \\%% & %8.0f \\%% & %8.1f \\%% & %8.0f \\%% \\\\ \n', 100*FullResults(1,1).Output.StationaryEqmStats(1).SSdivy, 100*FullResults(1,1).Output.StationaryEqmStats(2).SSdivy, 100*FullResults(2,1).Output.StationaryEqmStats(1).SSdivy, 100*FullResults(2,1).Output.StationaryEqmStats(2).SSdivy, 100*FullResults(3,1).Output.StationaryEqmStats(1).SSdivy, 100*FullResults(3,1).Output.StationaryEqmStats(2).SSdivy);
fprintf(FID, 'cv(lab)   & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(1,1).Output.StationaryEqmStats(1).CoeffOfVariance_l, FullResults(1,1).Output.StationaryEqmStats(2).CoeffOfVariance_l, FullResults(2,1).Output.StationaryEqmStats(1).CoeffOfVariance_l, FullResults(2,1).Output.StationaryEqmStats(2).CoeffOfVariance_l, FullResults(3,1).Output.StationaryEqmStats(1).CoeffOfVariance_l, FullResults(3,1).Output.StationaryEqmStats(2).CoeffOfVariance_l);
fprintf(FID, 'cv(weal)  & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(1,1).Output.StationaryEqmStats(1).CoeffOfVariance_a, FullResults(1,1).Output.StationaryEqmStats(2).CoeffOfVariance_a, FullResults(2,1).Output.StationaryEqmStats(1).CoeffOfVariance_a, FullResults(2,1).Output.StationaryEqmStats(2).CoeffOfVariance_a, FullResults(3,1).Output.StationaryEqmStats(1).CoeffOfVariance_a, FullResults(3,1).Output.StationaryEqmStats(2).CoeffOfVariance_a);
fprintf(FID, '$EV^{ss}$ & -- & %8.2f \\%% & -- & %8.2f \\%% & -- & %8.2f \\%% \\\\ \n',    100*FullResults(1,1).Output.StationaryEqmStats(2).EV_SS, 100*FullResults(2,1).Output.StationaryEqmStats(2).EV_SS, 100*FullResults(3,1).Output.StationaryEqmStats(2).EV_SS );
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: Do not attempt to replicate Co-worker mean as I do not know definition. Avg size of entering firms is defined in terms of nprime, while avg size of exiting firms is defined in terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Figures 1 & 2


figure(1)
age_axis_Jr=19+(1:1:Params.Jr); % Just working age
subplot(3,1,1); plot(age_axis_Jr, FullResults(1,1).Output.Figures_LifeCycleProfiles.HoursWorked_init(1:Params.Jr)) 
hold on
plot(age_axis_Jr, FullResults(1,1).Output.Figures_LifeCycleProfiles.HoursWorked_final(1:Params.Jr),'--')
hold off
title('No heterogeneity')
ylim([0,0.5])
subplot(3,1,2); plot(age_axis_Jr, FullResults(1,2).Output.Figures_LifeCycleProfiles.HoursWorked_init(1:Params.Jr)) 
hold on
plot(age_axis_Jr, FullResults(1,2).Output.Figures_LifeCycleProfiles.HoursWorked_final(1:Params.Jr),'--')
hold off
title('Symmetric Heterogeneity')
ylim([0,0.5])
subplot(3,1,3); plot(age_axis_Jr, FullResults(1,3).Output.Figures_LifeCycleProfiles.HoursWorked_init(1:Params.Jr)) 
hold on
plot(age_axis_Jr, FullResults(1,3).Output.Figures_LifeCycleProfiles.HoursWorked_final(1:Params.Jr),'--')
hold off
title('Asymmetric heterogeneity')
ylim([0,0.5])
xlabel('Age')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure1.png')


figure(2)
age_axis_J=19+(1:1:Params.J); % All ages
subplot(3,1,1); plot(age_axis_J, FullResults(1,1).Output.Figures_LifeCycleProfiles.Assets_init) 
hold on
plot(age_axis_J, FullResults(1,1).Output.Figures_LifeCycleProfiles.Assets_final,'--')
hold off
title('No heterogeneity')
ylim([0,15])
subplot(3,1,2); plot(age_axis_J, FullResults(1,2).Output.Figures_LifeCycleProfiles.Assets_init) 
hold on
plot(age_axis_J, FullResults(1,2).Output.Figures_LifeCycleProfiles.Assets_final,'--')
hold off
title('Symmetric Heterogeneity')
ylim([0,15])
subplot(3,1,3); plot(age_axis_J, FullResults(1,3).Output.Figures_LifeCycleProfiles.Assets_init) 
hold on
plot(age_axis_J, FullResults(1,3).Output.Figures_LifeCycleProfiles.Assets_final,'--')
hold off
title('Asymmetric heterogeneity')
ylim([0,15])
xlabel('Age')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure2.png')

%% Figures 3, 9, and 12
% Create graphs of output per capita, interest rates, capital-output ratio,
% and labor supply over the transition path.

figure(3)
subplot(2,2,1)
plot(0:1:T, FullResults(1,1).Output.Figures_TransitionAggregates.outputpercapita)
hold on
plot(0:1:T, FullResults(2,1).Output.Figures_TransitionAggregates.outputpercapita,'--')
plot(0:1:T, FullResults(3,1).Output.Figures_TransitionAggregates.outputpercapita,'-.')
hold off
title('Evolution of Output per capita')
ylim([1,2])
xlabel('Time')
subplot(2,2,2)
plot(0:1:T, FullResults(1,1).Output.Figures_TransitionAggregates.interestrate)
hold on
plot(0:1:T, FullResults(2,1).Output.Figures_TransitionAggregates.interestrate,'--')
plot(0:1:T, FullResults(3,1).Output.Figures_TransitionAggregates.interestrate,'-.')
hold off
title('Evolution of Interest Rate')
ylim([0,0.08])
xlabel('Time')
subplot(2,2,3)
plot(0:1:T, FullResults(1,1).Output.Figures_TransitionAggregates.capitaloutputratio)
hold on
plot(0:1:T, FullResults(2,1).Output.Figures_TransitionAggregates.capitaloutputratio,'--')
plot(0:1:T, FullResults(3,1).Output.Figures_TransitionAggregates.capitaloutputratio,'-.')
hold off
title('Evolution of Capital-Output Ratio')
ylim([2.5,5])
xlabel('Time')
subplot(2,2,4)
plot(0:1:T, 100*FullResults(1,1).Output.Figures_TransitionAggregates.hoursworked,'DisplayName','No Het')
hold on
plot(0:1:T, 100*FullResults(2,1).Output.Figures_TransitionAggregates.hoursworked,'--','DisplayName','Sym Het')
plot(0:1:T, 100*FullResults(3,1).Output.Figures_TransitionAggregates.hoursworked,'-.','DisplayName','Asym Het') % .hoursworked appears to be what CK1999, not .laborsupply
hold off
title('Evolution of Labor Supply (hours worked)')
% ylim([22,32])
xlabel('Time')
% Create a label for this fourth subplot but located it so it appears as if
% it is a superlegend: https://au.mathworks.com/matlabcentral/answers/387391-multiple-plots-with-same-legend#answer_309608
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.43;
Lgnd.Position(2) = 0.49;
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure3.png')


figure(9)
subplot(2,2,1)
plot(0:1:T, FullResults(1,2).Output.Figures_TransitionAggregates.outputpercapita)
hold on
plot(0:1:T, FullResults(2,2).Output.Figures_TransitionAggregates.outputpercapita,'--')
plot(0:1:T, FullResults(3,2).Output.Figures_TransitionAggregates.outputpercapita,'-.')
hold off
title('Evolution of Output per capita')
ylim([1,2])
xlabel('Time')
subplot(2,2,2)
plot(0:1:T, FullResults(1,2).Output.Figures_TransitionAggregates.interestrate)
hold on
plot(0:1:T, FullResults(2,2).Output.Figures_TransitionAggregates.interestrate,'--')
plot(0:1:T, FullResults(3,2).Output.Figures_TransitionAggregates.interestrate,'-.')
hold off
title('Evolution of Interest Rate')
ylim([0,0.08])
xlabel('Time')
subplot(2,2,3)
plot(0:1:T, FullResults(1,2).Output.Figures_TransitionAggregates.capitaloutputratio)
hold on
plot(0:1:T, FullResults(2,2).Output.Figures_TransitionAggregates.capitaloutputratio,'--')
plot(0:1:T, FullResults(3,2).Output.Figures_TransitionAggregates.capitaloutputratio,'-.')
hold off
title('Evolution of Capital-Output Ratio')
ylim([2.5,5])
xlabel('Time')
subplot(2,2,4)
plot(0:1:T, 100*FullResults(1,2).Output.Figures_TransitionAggregates.hoursworked,'DisplayName','No Het')
hold on
plot(0:1:T, 100*FullResults(2,2).Output.Figures_TransitionAggregates.hoursworked,'--','DisplayName','Sym Het')
plot(0:1:T, 100*FullResults(3,2).Output.Figures_TransitionAggregates.hoursworked,'-.','DisplayName','Asym Het') % .hoursworked appears to be what CK1999, not .laborsupply
hold off
title('Evolution of Labor Supply (hours worked)')
% ylim([22,32])
xlabel('Time')
% Create a label for this fourth subplot but located it so it appears as if
% it is a superlegend: https://au.mathworks.com/matlabcentral/answers/387391-multiple-plots-with-same-legend#answer_309608
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.43;
Lgnd.Position(2) = 0.49;
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure9.png')


figure(12)
subplot(2,2,1)
plot(0:1:T, FullResults(1,3).Output.Figures_TransitionAggregates.outputpercapita)
hold on
plot(0:1:T, FullResults(2,3).Output.Figures_TransitionAggregates.outputpercapita,'--')
plot(0:1:T, FullResults(3,3).Output.Figures_TransitionAggregates.outputpercapita,'-.')
hold off
title('Evolution of Output per capita')
ylim([1,2])
xlabel('Time')
subplot(2,2,2)
plot(0:1:T, FullResults(1,3).Output.Figures_TransitionAggregates.interestrate)
hold on
plot(0:1:T, FullResults(2,3).Output.Figures_TransitionAggregates.interestrate,'--')
plot(0:1:T, FullResults(3,3).Output.Figures_TransitionAggregates.interestrate,'-.')
hold off
title('Evolution of Interest Rate')
ylim([0,0.08])
xlabel('Time')
subplot(2,2,3)
plot(0:1:T, FullResults(1,3).Output.Figures_TransitionAggregates.capitaloutputratio)
hold on
plot(0:1:T, FullResults(2,3).Output.Figures_TransitionAggregates.capitaloutputratio,'--')
plot(0:1:T, FullResults(3,3).Output.Figures_TransitionAggregates.capitaloutputratio,'-.')
hold off
title('Evolution of Capital-Output Ratio')
ylim([2.5,5])
xlabel('Time')
subplot(2,2,4)
plot(0:1:T, 100*FullResults(1,3).Output.Figures_TransitionAggregates.hoursworked,'DisplayName','No Het')
hold on
plot(0:1:T, 100*FullResults(2,3).Output.Figures_TransitionAggregates.hoursworked,'--','DisplayName','Sym Het')
plot(0:1:T, 100*FullResults(3,3).Output.Figures_TransitionAggregates.hoursworked,'-.','DisplayName','Asym Het') % .hoursworked appears to be what CK1999, not .laborsupply
hold off
title('Evolution of Labor Supply (hours worked)')
% ylim([22,32])
xlabel('Time')
% Create a label for this fourth subplot but located it so it appears as if
% it is a superlegend: https://au.mathworks.com/matlabcentral/answers/387391-multiple-plots-with-same-legend#answer_309608
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.43;
Lgnd.Position(2) = 0.49;
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure12.png')

%% Figures 4, 10, 13: Vote shares
figure(4)
% Start with the 'population share' line.
populationshares=shiftdim(sum(sum(FullResults(1,1).Output.StationaryDist_init,1),2),2); % Note, in principle could just use mewj, but lets do this in following manner
plot(age_axis_J, populationshares,':');
% double-check: should equal Params.(AgeWeightsParamNames{1})
VoteInFavor_nohet=(FullResults(1,1).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
VoteInFavor_symhet=(FullResults(2,1).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
VoteInFavor_asymhet=(FullResults(3,1).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
ShareOfVotesInFavor_nohet=shiftdim(sum(sum(VoteInFavor_nohet.*FullResults(1,1).Output.StationaryDist_init,1),2),2);
ShareOfVotesInFavor_symhet=shiftdim(sum(sum(VoteInFavor_symhet.*FullResults(2,1).Output.StationaryDist_init,1),2),2);
ShareOfVotesInFavor_asymhet=shiftdim(sum(sum(VoteInFavor_asymhet.*FullResults(3,1).Output.StationaryDist_init,1),2),2);
hold on
plot(age_axis_J, ShareOfVotesInFavor_nohet, age_axis_J, ShareOfVotesInFavor_symhet,'--',age_axis_J, ShareOfVotesInFavor_asymhet,'-.')
hold off
legend('Pop. Share','No Het','Sym Het','Asym Het')
xlabel('age')
ylabel('Percentage of total votes')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure4.png')


figure(10)
% Start with the 'population share' line.
populationshares=shiftdim(sum(sum(FullResults(1,1).Output.StationaryDist_init,1),2),2);
plot(age_axis_J, populationshares,':');
% double-check: should equal Params.(AgeWeightsParamNames{1})
VoteInFavor_nohet=(FullResults(1,2).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
VoteInFavor_symhet=(FullResults(2,2).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
VoteInFavor_asymhet=(FullResults(3,2).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
ShareOfVotesInFavor_nohet=shiftdim(sum(sum(VoteInFavor_nohet.*FullResults(1,1).Output.StationaryDist_init,1),2),2);
ShareOfVotesInFavor_symhet=shiftdim(sum(sum(VoteInFavor_symhet.*FullResults(2,1).Output.StationaryDist_init,1),2),2);
ShareOfVotesInFavor_asymhet=shiftdim(sum(sum(VoteInFavor_asymhet.*FullResults(3,1).Output.StationaryDist_init,1),2),2);
hold on
plot(age_axis_J, ShareOfVotesInFavor_nohet, age_axis_J, ShareOfVotesInFavor_symhet,'--',age_axis_J, ShareOfVotesInFavor_asymhet,'-.')
hold off
legend('Pop. Share','No Het','Sym Het','Asym Het')
xlabel('age')
ylabel('Percentage of total votes')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure10.png')


figure(13)
% Start with the 'population share' line.
populationshares=shiftdim(sum(sum(FullResults(1,1).Output.StationaryDist_init,1),2),2);
plot(age_axis_J, populationshares,':');
% double-check: should equal Params.(AgeWeightsParamNames{1})
VoteInFavor_nohet=(FullResults(1,3).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
VoteInFavor_symhet=(FullResults(2,3).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
VoteInFavor_asymhet=(FullResults(3,3).Output.EquivVariation>=0); % agents who are indifferent are assumed to vote in favour
ShareOfVotesInFavor_nohet=shiftdim(sum(sum(VoteInFavor_nohet.*FullResults(1,1).Output.StationaryDist_init,1),2),2);
ShareOfVotesInFavor_symhet=shiftdim(sum(sum(VoteInFavor_symhet.*FullResults(2,1).Output.StationaryDist_init,1),2),2);
ShareOfVotesInFavor_asymhet=shiftdim(sum(sum(VoteInFavor_asymhet.*FullResults(3,1).Output.StationaryDist_init,1),2),2);
hold on
plot(age_axis_J, ShareOfVotesInFavor_nohet, age_axis_J, ShareOfVotesInFavor_symhet,'--',age_axis_J, ShareOfVotesInFavor_asymhet,'-.')
hold off
legend('Pop. Share','No Het','Sym Het','Asym Het')
xlabel('age')
ylabel('Percentage of total votes')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure13.png')


%% Figures 5,6,7

figure(5)
EV=FullResults(1,1).Output.EquivVariation;
plot(a_grid, EV(:,1,20-19),'r', a_grid, EV(:,1,30-19),'b', a_grid, EV(:,1,60-19),'g', a_grid, zeros(size(a_grid)),'k')
legend('Age 20', 'Age 30', 'Age 60')
title('Welfare effects of Reform A: No heterogeneity')
ylim([-0.5,0.2])
xlabel('Asset Position')
ylabel('Cons. Equiv. Var.')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure5.png')


figure(6)
EV=FullResults(2,1).Output.EquivVariation;
hold on
Y1=plot(a_grid, EV(:,1,20-19),'r');
Y2=plot(a_grid, EV(:,1,30-19),'b');
Y3=plot(a_grid, EV(:,1,60-19),'g');
plot( a_grid, zeros(size(a_grid)),'k');
plot(a_grid, EV(:,2,20-19),'--r', a_grid, EV(:,2,30-19),'--b', a_grid, EV(:,2,60-19),'--g')
% Requires setting up hidden lines with the appropriate style
H1 = plot(a_grid,EV(:,1,20-19), '-', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % solid line (invisible and black)
H2 = plot(a_grid,EV(:,1,20-19), '--', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % dashed line (invisible and black)
hold off
legend([Y1,Y2,Y3,H1,H2],'Age 20', 'Age 30', 'Age 60','Bad Shock','Good Shock');
title('Welfare effects of Reform A: Symm. heterogeneity')
ylim([-0.5,0.2])
xlabel('Asset Position')
ylabel('Cons. Equiv. Var.')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure6.png')


figure(7)
EV=FullResults(3,1).Output.EquivVariation;
hold on
Y1=plot(a_grid, EV(:,1,20-19),'r');
Y2=plot(a_grid, EV(:,1,30-19),'b');
Y3=plot(a_grid, EV(:,1,60-19),'g');
plot(a_grid, zeros(size(a_grid)),'k');
plot(a_grid, EV(:,2,20-19),'--r', a_grid, EV(:,2,30-19),'--b', a_grid, EV(:,2,60-19),'--g')
% Requires setting up hidden lines with the appropriate style
H1 = plot(a_grid,EV(:,1,20-19), '-', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % solid line (invisible and black)
H2 = plot(a_grid,EV(:,1,20-19), '--', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % dashed line (invisible and black)
hold off
legend([Y1,Y2,Y3,H1,H2],'Age 20', 'Age 30', 'Age 60','Bad Shock','Good Shock');
title('Welfare effects of Reform A: Asymm. heterogeneity')
ylim([-0.8,0.2])
xlabel('Asset Position')
ylabel('Cons. Equiv. Var.')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure7.png')



%% Figure 8

figure(8)
EV=FullResults(3,1).Output.EquivVariation_FixedPrices;
hold on
Y1=plot(a_grid, EV(:,1,20-19),'r');
Y2=plot(a_grid, EV(:,1,30-19),'b');
Y3=plot(a_grid, EV(:,1,60-19),'g');
plot(a_grid, zeros(size(a_grid)),'k')
plot(a_grid, EV(:,2,20-19),'--r', a_grid, EV(:,2,30-19),'--b', a_grid, EV(:,2,60-19),'--g')
% Requires setting up hidden lines with the appropriate style
H1 = plot(a_grid,EV(:,1,20-19), '-', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % solid line (invisible and black)
H2 = plot(a_grid,EV(:,1,20-19), '--', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % dashed line (invisible and black)
hold off
legend([Y1,Y2,Y3,H1,H2],'Age 20', 'Age 30', 'Age 60','Bad Shock','Good Shock');
title('Welfare effects of Reform A: Asymm. heterogeneity')
ylim([-0.8,0.2])
xlabel('Asset Position')
ylabel('Cons. Equiv. Var.')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure8.png')


%% Figure 11

% NOTE: Uses different ages to all other figures!
figure(11)
EV=FullResults(3,2).Output.EquivVariation;
hold on
Y1=plot(a_grid, EV(:,1,20-19),'r');
Y2=plot(a_grid, EV(:,1,45-19),'b');
Y3=plot(a_grid, EV(:,1,81-19),'g');
plot(a_grid, zeros(size(a_grid)),'k')
plot(a_grid, EV(:,2,20-19),'--r', a_grid, EV(:,2,45-19),'--b', a_grid, EV(:,2,81-19),'--g')
% Requires setting up hidden lines with the appropriate style
H1 = plot(a_grid,EV(:,1,20-19), '-', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % solid line (invisible and black)
H2 = plot(a_grid,EV(:,1,20-19), '--', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % dashed line (invisible and black)
hold off
legend([Y1,Y2,Y3,H1,H2],'Age 20', 'Age 45', 'Age 81','Bad Shock','Good Shock');
title('Welfare effects of Reform B: Asymm. heterogeneity')
xlabel('Asset Position')
ylabel('Cons. Equiv. Var.')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure11.png')


%% Table 6

% Reform A
% VoteShareInFavor_sym2state=sum(sum(sum((FullResults(2,1).Output.EquivVariation>=0).*FullResults(2,1).Output.StationaryDist_init))); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state=sum(sum(sum((FullResults(4,1).Output.EquivVariation>=0).*FullResults(2,1).Output.StationaryDist_init))); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym3state=sum(sum(sum((FullResults(5,1).Output.EquivVariation>=0).*FullResults(5,1).Output.StationaryDist_init))); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym5state=sum(sum(sum((FullResults(6,1).Output.EquivVariation>=0).*FullResults(6,1).Output.StationaryDist_init))); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_asym2state=sum(sum(sum((FullResults(3,1).Output.EquivVariation>=0).*FullResults(3,1).Output.StationaryDist_init))); % agents who are indifferent are assumed to vote in favour

%Table 6
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table6.tex', 'w');
fprintf(FID, 'Dispersion of labor earnings, wealth, votes for Reform A \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \\hline \n');
fprintf(FID, ' & Sym. 2 States &  Sym. 3 States &  Sym. 5 States &  Asym. 2 States   \\\\ \\hline \n');
fprintf(FID, 'cv(lab)   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(4,1).Output.StationaryEqmStats(1).CoeffOfVariance_l, FullResults(5,1).Output.StationaryEqmStats(2).CoeffOfVariance_l, FullResults(6,1).Output.StationaryEqmStats(1).CoeffOfVariance_l, FullResults(3,1).Output.StationaryEqmStats(2).CoeffOfVariance_l);
fprintf(FID, 'cv(weal)  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(4,1).Output.StationaryEqmStats(1).CoeffOfVariance_a, FullResults(5,1).Output.StationaryEqmStats(2).CoeffOfVariance_a, FullResults(6,1).Output.StationaryEqmStats(1).CoeffOfVariance_a, FullResults(3,1).Output.StationaryEqmStats(2).CoeffOfVariance_a);
fprintf(FID, 'Votes  & %8.2f \\%% & %8.2f  \\%% & %8.2f  \\%% & %8.2f  \\%% \\\\ \n', 100*VoteShareInFavor_sym2state, 100*VoteShareInFavor_sym3state, 100*VoteShareInFavor_sym5state, 100*VoteShareInFavor_asym2state);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: notice. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Table 7

% Reform A
VoteShareInFavor_sym2state_type1=sum(sum(sum((FullResults(7,1).Output.EquivVariation.type1>=0).*FullResults(7,1).Output.StationaryDist_init.type1)))*FullResults(7,1).Output.StationaryDist_init.ptweights(1); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state=VoteShareInFavor_sym2state_type1;
VoteShareInFavor_sym2state2type_type1=sum(sum(sum(sum((FullResults(8,1).Output.EquivVariation.type1>=0).*FullResults(8,1).Output.StationaryDist_init.type1))))*FullResults(8,1).Output.StationaryDist_init.ptweights(1); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state2type_type2=sum(sum(sum(sum((FullResults(8,1).Output.EquivVariation.type2>=0).*FullResults(8,1).Output.StationaryDist_init.type2))))*FullResults(8,1).Output.StationaryDist_init.ptweights(2); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state2type=VoteShareInFavor_sym2state2type_type1+VoteShareInFavor_sym2state2type_type2; % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state3type_type1=sum(sum(sum(sum((FullResults(9,1).Output.EquivVariation.type1>=0).*FullResults(9,1).Output.StationaryDist_init.type1))))*FullResults(9,1).Output.StationaryDist_init.ptweights(1); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state3type_type2=sum(sum(sum(sum((FullResults(9,1).Output.EquivVariation.type2>=0).*FullResults(9,1).Output.StationaryDist_init.type2))))*FullResults(9,1).Output.StationaryDist_init.ptweights(2); % agents who are indifferent are assumed to vote in favour
VoteShareInFavor_sym2state3type=VoteShareInFavor_sym2state3type_type1+VoteShareInFavor_sym2state3type_type2; % agents who are indifferent are assumed to vote in favour

%Table 7
FID = fopen('./SavedOutput/LatexInputs/ConesaKrueger1999_Table7.tex', 'w');
fprintf(FID, 'Dispersion of labor earnings, wealth, votes for Reform A \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccc} \n \\hline \\hline \n');
fprintf(FID, ' & Sym. 2 States  & Sym. 2 type, 2 States  & Sym. 3 types, 2 States  \\\\ \\hline \n');
fprintf(FID, 'cv(lab)   & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(7,1).Output.StationaryEqmStats(1).CoeffOfVariance_l, FullResults(8,1).Output.StationaryEqmStats(2).CoeffOfVariance_l, FullResults(9,1).Output.StationaryEqmStats(1).CoeffOfVariance_l);
fprintf(FID, 'cv(weal)  & %8.2f & %8.2f & %8.2f \\\\ \n',    FullResults(7,1).Output.StationaryEqmStats(1).CoeffOfVariance_a, FullResults(8,1).Output.StationaryEqmStats(2).CoeffOfVariance_a, FullResults(9,1).Output.StationaryEqmStats(1).CoeffOfVariance_a);
fprintf(FID, 'Votes  & %8.2f  \\%% & %8.2f  \\%% & %8.2f  \\%% \\\\ \n', 100*VoteShareInFavor_sym2state, 100*VoteShareInFavor_sym2state2type, 100*VoteShareInFavor_sym2state3type);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: terms of n. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Figure 14

% Symmetric, 2 states, 3 types, reform A.

% Unclear what to do with age dimension? I tried taking means of EV over the age
% dimension to eliminate it but this gave weird answers. It is unclear what CK1999 did. Does not appear to
% be any mention in text (their previous similar Figures are all for
% specific ages). I have therefore ended up just drawing the age 20
% version as out of the age 20, 30 & 60 versions this was the one that
% looked most like that of CK1999.

clear EV
figure(14)
EV.type1=FullResults(9,1).Output.EquivVariation.type1;
EV.type2=FullResults(9,1).Output.EquivVariation.type2;
EV.type3=FullResults(9,1).Output.EquivVariation.type3;
% EV=sum(EV.*FullResults(9,1).Output.StationaryDist_init, 4)/N_j; % Take mean over age distriubtion
EV.type1=EV.type1(:,:,:,20-19); % Alternative, just pick a specific age
EV.type2=EV.type2(:,:,:,20-19); % Alternative, just pick a specific age
EV.type3=EV.type3(:,:,:,20-19); % Alternative, just pick a specific age
hold on
Y1=plot(a_grid, EV.type1(:,1),'r');
Y2=plot(a_grid, EV.type2(:,1),'b');
Y3=plot(a_grid, EV.type3(:,1),'g');
plot(a_grid, zeros(size(a_grid)),'k') % Bad shock
plot(a_grid, EV.type1(:,2),'--r', a_grid, EV.type2(:,2),'--b', a_grid, EV.type3(:,2),'--g') % Good shock
% Requires setting up hidden lines with the appropriate style
H1 = plot(a_grid,EV.type1(:,1), '-', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % solid line (invisible and black)
H2 = plot(a_grid,EV.type1(:,1), '--', 'LineWidth', 2, 'Color', 'k', 'Visible', 'off'); % dashed line (invisible and black)
hold off
legend([Y1,Y2,Y3,H1,H2],'Type 1 (Age 20)', 'Type 2 (Age 20)', 'Type 3 (Age 20)','Bad Shock','Good Shock');
title('Welfare effects of Reform A with permanent shocks')
xlabel('Asset Position')
ylabel('Cons. Equiv. Var.')
saveas(gcf,'./SavedOutput/Graphs/ConesaKrueger1999_Figure14.png')










