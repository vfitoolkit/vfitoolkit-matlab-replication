%% Replication of Imrohoroglu, Imrohoroglu & Joines (1995) - A life cycle analysis of social security

% A few lines I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

% No decision variable: IIJ1995 has inelastic labor supply
% One endogenous state variable: assets (k)
% One stochastic exogenous state variable (z, which is either 'employed' or 'unemployed')
% Age

Params.J=65; % Number of period in life-cycle (represents ages 21 to 85)

% Grid sizes to use
% % labour supply is exogenous (IIJ1995 refers to this as 'perfectly inelastic labor supply'
n_k=1251; % Assets
n_z=2; % employed or unemployed (IIJ1995 refers to this shock as Markow, but in calibration actually imposes that it is iid).

%% Replication notes:
% Figures 6 & 7, IIJ1995 do not explicitly state which social security
% replacement rate is used for these. Since 0.3 is often used for other
% Figures. Also unclear if Figure 7 plots age-conditonal pdfs, or if it
% incorporates the different mass of agents of each age. I plot the age
% conditional pdfs.

%% The values to be looped over

beta_vec=[0.98,1.011]; % Discount factor
b_vec=0:0.1:1; % Social security benefit rate
n_vec=[0,0.012]; % Population growth rate
gamma_vec=[1.5,2,4];
g_vec=[0,0.022];
% stochprobofdeath (what IIJ1995 call certain and uncertain lifetimes)
% MedicalShock

% beta_c=1; b_c=11; n_c=1; gamma_c=3; g_c=1; MedicalShock=1; stochprobofdeath=0;

% I save the results in two halves, just in case there are problems halfway through.
% Following two (commented out) lines allow me to resume
% load ./SavedOutput/IIJ1995_FullOutput1.mat FullOutput1
load ./SavedOutput/IIJ1995_FullOutput2.mat FullOutput2

for beta_c=2:2 % 1:2
    Params.beta=beta_vec(beta_c);
    for b_c=1:11 % 7
        Params.b=b_vec(b_c);
        for n_c=2:2 %1:2
            Params.n=n_vec(n_c);
            for gamma_c=2:2 % 1:3
                Params.gamma=gamma_vec(gamma_c);
                for g_c=1:1 % 1:2
                    Params.g=g_vec(g_c);
                    for MedicalShock=0:0 % 0:2
                        Params.MedicalShock=MedicalShock;
                        for stochprobofdeath=1:1 %0:1
                            fprintf('Currently doing: %i %i %i %i %i %i %i \n ', beta_c, b_c, n_c, gamma_c, g_c, MedicalShock, stochprobofdeath)
                            if beta_c==2 && b_c>1 && n_c==2 && gamma_c==2 && g_c==1 && MedicalShock==0 && stochprobofdeath==1
                                calcwelfarebenefits=1; 
                                % Only calculate the welfare benefits when needed (for Table 3).
                                % Calculating welfare benefits takes much more time than everything else, hence only do it when necessary.
                            else
                                calcwelfarebenefits=0;
                            end
                            
                            Output=ImrohorogluImrohorogluJoines1995_Fn(Params,n_k,n_z,stochprobofdeath,calcwelfarebenefits);
                            if beta_c==1
                                FullOutput1(beta_c,b_c,n_c,gamma_c,g_c,MedicalShock+1,stochprobofdeath+1).Output=Output;
                                save ./SavedOutput/IIJ1995_FullOutput1.mat FullOutput1 -v7.3
                            elseif beta_c==2
                                FullOutput2(beta_c,b_c,n_c,gamma_c,g_c,MedicalShock+1,stochprobofdeath+1).Output=Output;
                                save ./SavedOutput/IIJ1995_FullOutput2.mat FullOutput2 -v7.3
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Create Figures

load ./SavedOutput/IIJ1995_FullOutput1.mat FullOutput1
load ./SavedOutput/IIJ1995_FullOutput2.mat FullOutput2

% Figure 1
figure(1)
Params=FullOutput1(1,1,1,1,1,1,1).Output.Params;
discountfactor=cumprod(((Params.beta).*ones(1,Params.J)).*Params.sj);
subplot(2,2,1); plot(1:1:Params.J,discountfactor)
ylim([0,2.4])
title('beta=0.98, Certain Lifetimes')
Params=FullOutput1(1,1,1,1,1,1,2).Output.Params;
discountfactor=cumprod(((Params.beta).*ones(1,Params.J)).*Params.sj);
subplot(2,2,2); plot(1:1:Params.J,discountfactor)
ylim([0,2.4])
title('beta=0.98, Uncertain Lifetimes')
Params=FullOutput2(2,1,1,1,1,1,1).Output.Params;
discountfactor=cumprod(((Params.beta).*ones(1,Params.J)).*Params.sj);
subplot(2,2,3); plot(1:1:Params.J,discountfactor)
ylim([0,2.4])
title('beta=1.011, Certain Lifetimes')
Params=FullOutput2(2,1,1,1,1,1,2).Output.Params;
discountfactor=cumprod(((Params.beta).*ones(1,Params.J)).*Params.sj);
subplot(2,2,4); plot(1:1:Params.J,discountfactor)
ylim([0,2.4])
title('beta=1.011, Uncertain Lifetimes')
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure1.png')

% Figure 2
Params=FullOutput2(2,1,2,2,1,1,2).Output.Params;
figure(2)
plot(1:1:Params.J,100*FullOutput2(2,1,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean/FullOutput2(2,1,2,2,1,1,2).Output.C) % The 3 is consumption
hold on
plot(1:1:Params.J,100*FullOutput2(2,4,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean/FullOutput2(2,4,2,2,1,1,2).Output.C)
plot(1:1:Params.J,100*FullOutput2(2,7,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean/FullOutput2(2,7,2,2,1,1,2).Output.C)
% The 'planners' can be calculated directly:
% IIJ1995 bottom pg 96, top pg 97
% Consumption Euler eqn: c_t^-gamma=(1+r)*beta*c_{t+1}^-gamma. This defines the 'planners' life-cycle profile for consumption (which provides perfect insurance against employment and mortality risks).
% Rearranging gives (c_{t+1}/c_t)=(1/((1+r)*beta))^(-1/gamma)
% So now we just need to find the planner interest rate and we will have the answer
% IIJ1995 top pg 97 explain how this is just the same exercise as finding
% the 'golden rule' level of capital in the Solow growth model.
% Look at any notes on 'golden rule' in Solow model and it will explain how to derive: KdivN_G=((n+g+delta)/alpha)^(-1/(1-alpha));
KdivN_G=((Params.n+Params.g+Params.delta)/Params.alpha)^(-1/(1-Params.alpha));
% We can then get the planners interest rate from r+delta=MPK (marginal product of capital)
planner_interestrate=Params.alpha*Params.A*(KdivN_G^(Params.alpha-1))-Params.delta;
consumptiongrowthrate=cumprod(((1/((1+planner_interestrate)*(FullOutput1(1,1,2,2,1,1,2).Output.Params.beta)))^(-1/Params.gamma))*ones(Params.J,1));
plannersconsumptionprofile=100*consumptiongrowthrate/sum(FullOutput2(2,1,2,2,1,1,2).Output.Params.mewj.*consumptiongrowthrate'); % Normalize by aggregate consumption.
plot(1:1:Params.J, plannersconsumptionprofile) % Consumption Euler eqn: c_t^-gamma=(1+r)*beta*c_{t+1}^-gamma. This defines the 'planners' life-cycle profile for consumption (which provides perfect insurance against employment and mortality risks).
legend('b=0', 'b=0.3', 'b=0.6', 'planner')
hold off
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure2.png')

% Figure 3
figure(3)
% Social security replacement rate of 30%
plot(1:1:Params.J,FullOutput2(2,4,2,2,1,1,2).Output.LifeCycleProfiles(4).Mean) % The 4 is income
hold on
plot(1:1:Params.J,FullOutput2(2,4,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean) % The 3 is consumption
legend('income','consumption')
hold off
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure3.png')

% Figure 4
figure(4)
% The 'CES data' is not plotted as that is dependent on empirical data
plot(1:1:Params.J,FullOutput2(2,1,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean) % The 3 is consumption
hold on
plot(1:1:Params.J,FullOutput2(2,4,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean) % The 3 is consumption
plot(1:1:Params.J,FullOutput2(2,7,2,2,1,1,2).Output.LifeCycleProfiles(3).Mean) % The 3 is consumption
legend('b=0', 'b=0.3', 'b=0.6')
hold off
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure4.png')

% Figure 5
figure(5)
plot(1:1:Params.J,FullOutput2(2,1,2,2,1,1,2).Output.LifeCycleProfiles(1).Mean) % The 1 is assets
hold on
plot(1:1:Params.J,FullOutput2(2,4,2,2,1,1,2).Output.LifeCycleProfiles(1).Mean) % The 1 is assets
plot(1:1:Params.J,FullOutput2(2,7,2,2,1,1,2).Output.LifeCycleProfiles(1).Mean) % The 1 is assets
legend('b=0', 'b=0.3', 'b=0.6')
hold off
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure5.png')


% Figures 6 & 7, IIJ1995 do not explicitly state which social security
% replacement rate is used for these. Since 0.3 is often used for other
% Figures. Also unclear if Figure 7 plots age-conditonal pdfs, or if it
% incorporates the different mass of agents of each age. I plot the age
% conditional pdfs.

% Figure 6
figure(6)
AssetValuesOnGrid=shiftdim(FullOutput2(2,1,2,2,1,1,2).Output.ValuesOnGrid(1,:,:,:),1); % ValuesOnGrid(1,a,z,j), 1 is assets
StationaryDist=FullOutput2(2,1,2,2,1,1,2).Output.StationaryDist;
FractionOfHHsPerAssetPartition=zeros(2*15,1); % 15 is the max value of the asset grid
for ii=1:(2*15) % 15 is the max value of the asset grid
    Partition=logical((AssetValuesOnGrid>=((ii-1)/2)).*(AssetValuesOnGrid<((ii-1)/2+0.5)));
    FractionOfHHsPerAssetPartition(ii)=sum(StationaryDist(Partition));
end
plot(0.25:0.5:(15-0.25),FractionOfHHsPerAssetPartition)
xticks(0:1:15)
title('Assets')
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure6.png')

% a_grid only depends on n_k so is the same for everything
a_grid=FullOutput1(1,1,1,1,1,1,1).Output.a_grid;

% Figure 7
figure(7)
plot(a_grid,sum(StationaryDist(:,:,60),2)/sum(sum(StationaryDist(:,:,60))))
hold on
plot(a_grid,sum(StationaryDist(:,:,54),2)/sum(sum(StationaryDist(:,:,54))))
plot(a_grid,sum(StationaryDist(:,:,50),2)/sum(sum(StationaryDist(:,:,50))))
plot(a_grid,sum(StationaryDist(:,:,44),2)/sum(sum(StationaryDist(:,:,44))))
hold off
xticks(0:1:15)
legend('age 60','age 54','age 50','age 44')
title('Age-Conditional pdf of Assets')
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure7.png')


% Figure 7, alternative version using cdf (pdf are sensitive to grid sizes
% and so typically not a good thing to plot, unless you use a kernel pdf
% estimator to fit and plot them.)
figure(7)
plot(a_grid,cumsum(sum(StationaryDist(:,:,60),2))/sum(sum(StationaryDist(:,:,60))))
hold on
plot(a_grid,cumsum(sum(StationaryDist(:,:,54),2))/sum(sum(StationaryDist(:,:,54))))
plot(a_grid,cumsum(sum(StationaryDist(:,:,50),2))/sum(sum(StationaryDist(:,:,50))))
plot(a_grid,cumsum(sum(StationaryDist(:,:,44),2))/sum(sum(StationaryDist(:,:,44))))
hold off
xticks(0:1:15)
legend('age 60','age 54','age 50','age 44')
title('Age-Conditional cdf of Assets')
saveas(gcf,'./SavedOutput/Graphs/IIJ1995_Figure7cdf.png')


%% Create Tables

% Table 1
FID = fopen('./SavedOutput/LatexInputs/ImrohorogluImrohorogluJoines1995_Table1.tex', 'w');
fprintf(FID, 'Population growth and lifetime uncertainty, $\\beta=1.011$, $\\gamma = 2$ \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllllll} \n \\hline \\hline \n');
fprintf(FID, ' & Tax & Wage & Return to & Average & Capital & Average & Average  \\\\ \\hline \n');
fprintf(FID, ' b & rate & rate & capital & consumption & stock & income & utility  \\\\ \\hline \n');
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,1,2,2,1,1,2).Output.Params.b,FullOutput2(2,1,2,2,1,1,2).Output.tau_s,FullOutput2(2,1,2,2,1,1,2).Output.w,FullOutput2(2,1,2,2,1,1,2).Output.r,FullOutput2(2,1,2,2,1,1,2).Output.C,FullOutput2(2,1,2,2,1,1,2).Output.K,FullOutput2(2,1,2,2,1,1,2).Output.Q,FullOutput2(2,1,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.Params.b,FullOutput2(2,2,2,2,1,1,2).Output.tau_s,FullOutput2(2,2,2,2,1,1,2).Output.w,FullOutput2(2,2,2,2,1,1,2).Output.r,FullOutput2(2,2,2,2,1,1,2).Output.C,FullOutput2(2,2,2,2,1,1,2).Output.K,FullOutput2(2,2,2,2,1,1,2).Output.Q,FullOutput2(2,2,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,3,2,2,1,1,2).Output.Params.b,FullOutput2(2,3,2,2,1,1,2).Output.tau_s,FullOutput2(2,3,2,2,1,1,2).Output.w,FullOutput2(2,3,2,2,1,1,2).Output.r,FullOutput2(2,3,2,2,1,1,2).Output.C,FullOutput2(2,3,2,2,1,1,2).Output.K,FullOutput2(2,3,2,2,1,1,2).Output.Q,FullOutput2(2,3,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,4,2,2,1,1,2).Output.Params.b,FullOutput2(2,4,2,2,1,1,2).Output.tau_s,FullOutput2(2,4,2,2,1,1,2).Output.w,FullOutput2(2,4,2,2,1,1,2).Output.r,FullOutput2(2,4,2,2,1,1,2).Output.C,FullOutput2(2,4,2,2,1,1,2).Output.K,FullOutput2(2,4,2,2,1,1,2).Output.Q,FullOutput2(2,4,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,5,2,2,1,1,2).Output.Params.b,FullOutput2(2,5,2,2,1,1,2).Output.tau_s,FullOutput2(2,5,2,2,1,1,2).Output.w,FullOutput2(2,5,2,2,1,1,2).Output.r,FullOutput2(2,5,2,2,1,1,2).Output.C,FullOutput2(2,5,2,2,1,1,2).Output.K,FullOutput2(2,5,2,2,1,1,2).Output.Q,FullOutput2(2,5,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,6,2,2,1,1,2).Output.Params.b,FullOutput2(2,6,2,2,1,1,2).Output.tau_s,FullOutput2(2,6,2,2,1,1,2).Output.w,FullOutput2(2,6,2,2,1,1,2).Output.r,FullOutput2(2,6,2,2,1,1,2).Output.C,FullOutput2(2,6,2,2,1,1,2).Output.K,FullOutput2(2,6,2,2,1,1,2).Output.Q,FullOutput2(2,6,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,7,2,2,1,1,2).Output.Params.b,FullOutput2(2,7,2,2,1,1,2).Output.tau_s,FullOutput2(2,7,2,2,1,1,2).Output.w,FullOutput2(2,7,2,2,1,1,2).Output.r,FullOutput2(2,7,2,2,1,1,2).Output.C,FullOutput2(2,7,2,2,1,1,2).Output.K,FullOutput2(2,7,2,2,1,1,2).Output.Q,FullOutput2(2,7,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f  \\\\ \n',FullOutput2(2,11,2,2,1,1,2).Output.Params.b,FullOutput2(2,11,2,2,1,1,2).Output.tau_s,FullOutput2(2,11,2,2,1,1,2).Output.w,FullOutput2(2,11,2,2,1,1,2).Output.r,FullOutput2(2,11,2,2,1,1,2).Output.C,FullOutput2(2,11,2,2,1,1,2).Output.K,FullOutput2(2,11,2,2,1,1,2).Output.Q,FullOutput2(2,11,2,2,1,1,2).Output.Omega1);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 2
FID = fopen('./SavedOutput/LatexInputs/ImrohorogluImrohorogluJoines1995_Table2.tex', 'w');
fprintf(FID, 'The role of population growth and lifetime uncertainty \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllllllll} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{3}{l}{Zero population growth} & \\multicolumn{3}{l}{Population growth} & \\multicolumn{3}{l}{Population growth}  \\\\ \n');
fprintf(FID, ' & \\multicolumn{3}{l}{certain lifetimes} & \\multicolumn{3}{l}{certain lifetimes} & \\multicolumn{3}{l}{lifetime uncertainty}  \\\\  \\cline{2-4} \\cline{5-7} \\cline{8-10} \n');
fprintf(FID, ' b & K/Q & r & Utility & K/Q & r & Utility & K/Q & r & Utility  \\\\ \\hline \n');
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,1,2,2,1,1,2).Output.Params.b,FullOutput2(2,1,1,2,1,1,1).Output.KdivQ,FullOutput2(2,1,1,2,1,1,1).Output.r,FullOutput2(2,1,1,2,1,1,1).Output.Omega1,FullOutput2(2,1,2,2,1,1,1).Output.KdivQ,FullOutput2(2,1,2,2,1,1,1).Output.r,FullOutput2(2,1,2,2,1,1,1).Output.Omega1,FullOutput2(2,1,2,2,1,1,2).Output.KdivQ,FullOutput2(2,1,2,2,1,1,2).Output.r,FullOutput2(2,1,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.Params.b,FullOutput2(2,2,1,2,1,1,1).Output.KdivQ,FullOutput2(2,2,1,2,1,1,1).Output.r,FullOutput2(2,2,1,2,1,1,1).Output.Omega1,FullOutput2(2,2,2,2,1,1,1).Output.KdivQ,FullOutput2(2,2,2,2,1,1,1).Output.r,FullOutput2(2,2,2,2,1,1,1).Output.Omega1,FullOutput2(2,2,2,2,1,1,2).Output.KdivQ,FullOutput2(2,2,2,2,1,1,2).Output.r,FullOutput2(2,2,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,3,2,2,1,1,2).Output.Params.b,FullOutput2(2,3,1,2,1,1,1).Output.KdivQ,FullOutput2(2,3,1,2,1,1,1).Output.r,FullOutput2(2,3,1,2,1,1,1).Output.Omega1,FullOutput2(2,3,2,2,1,1,1).Output.KdivQ,FullOutput2(2,3,2,2,1,1,1).Output.r,FullOutput2(2,3,2,2,1,1,1).Output.Omega1,FullOutput2(2,3,2,2,1,1,2).Output.KdivQ,FullOutput2(2,3,2,2,1,1,2).Output.r,FullOutput2(2,3,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,4,2,2,1,1,2).Output.Params.b,FullOutput2(2,4,1,2,1,1,1).Output.KdivQ,FullOutput2(2,4,1,2,1,1,1).Output.r,FullOutput2(2,4,1,2,1,1,1).Output.Omega1,FullOutput2(2,4,2,2,1,1,1).Output.KdivQ,FullOutput2(2,4,2,2,1,1,1).Output.r,FullOutput2(2,4,2,2,1,1,1).Output.Omega1,FullOutput2(2,4,2,2,1,1,2).Output.KdivQ,FullOutput2(2,4,2,2,1,1,2).Output.r,FullOutput2(2,4,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,5,2,2,1,1,2).Output.Params.b,FullOutput2(2,5,1,2,1,1,1).Output.KdivQ,FullOutput2(2,5,1,2,1,1,1).Output.r,FullOutput2(2,5,1,2,1,1,1).Output.Omega1,FullOutput2(2,5,2,2,1,1,1).Output.KdivQ,FullOutput2(2,5,2,2,1,1,1).Output.r,FullOutput2(2,5,2,2,1,1,1).Output.Omega1,FullOutput2(2,5,2,2,1,1,2).Output.KdivQ,FullOutput2(2,5,2,2,1,1,2).Output.r,FullOutput2(2,5,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,6,2,2,1,1,2).Output.Params.b,FullOutput2(2,6,1,2,1,1,1).Output.KdivQ,FullOutput2(2,6,1,2,1,1,1).Output.r,FullOutput2(2,6,1,2,1,1,1).Output.Omega1,FullOutput2(2,6,2,2,1,1,1).Output.KdivQ,FullOutput2(2,6,2,2,1,1,1).Output.r,FullOutput2(2,6,2,2,1,1,1).Output.Omega1,FullOutput2(2,6,2,2,1,1,2).Output.KdivQ,FullOutput2(2,6,2,2,1,1,2).Output.r,FullOutput2(2,6,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,7,2,2,1,1,2).Output.Params.b,FullOutput2(2,7,1,2,1,1,1).Output.KdivQ,FullOutput2(2,7,1,2,1,1,1).Output.r,FullOutput2(2,7,1,2,1,1,1).Output.Omega1,FullOutput2(2,7,2,2,1,1,1).Output.KdivQ,FullOutput2(2,7,2,2,1,1,1).Output.r,FullOutput2(2,7,2,2,1,1,1).Output.Omega1,FullOutput2(2,7,2,2,1,1,2).Output.KdivQ,FullOutput2(2,7,2,2,1,1,2).Output.r,FullOutput2(2,7,2,2,1,1,2).Output.Omega1);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 3
FID = fopen('./SavedOutput/LatexInputs/ImrohorogluImrohorogluJoines1995_Table3.tex', 'w');
fprintf(FID, 'Welfare benefits of social security \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllllll} \n \\hline \\hline \n');
fprintf(FID, ' b         & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.1f \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.Params.b,FullOutput2(2,3,2,2,1,1,2).Output.Params.b,FullOutput2(2,4,2,2,1,1,2).Output.Params.b,FullOutput2(2,5,2,2,1,1,2).Output.Params.b,FullOutput2(2,6,2,2,1,1,2).Output.Params.b,FullOutput2(2,7,2,2,1,1,2).Output.Params.b,FullOutput2(2,11,2,2,1,1,2).Output.Params.b);
fprintf(FID, ' $\\kappa$ & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.kappa,FullOutput2(2,3,2,2,1,1,2).Output.kappa,FullOutput2(2,4,2,2,1,1,2).Output.kappa,FullOutput2(2,5,2,2,1,1,2).Output.kappa,FullOutput2(2,6,2,2,1,1,2).Output.kappa,FullOutput2(2,7,2,2,1,1,2).Output.kappa,FullOutput2(2,11,2,2,1,1,2).Output.kappa);
% fprintf(FID, ' $\\kappa$ & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\ \n',100*FullOutput2(2,2,2,2,1,1,2).Output.kappa,100*FullOutput2(2,3,2,2,1,1,2).Output.kappa,100*FullOutput2(2,4,2,2,1,1,2).Output.kappa,100*FullOutput2(2,5,2,2,1,1,2).Output.kappa,100*FullOutput2(2,6,2,2,1,1,2).Output.kappa,100*FullOutput2(2,7,2,2,1,1,2).Output.kappa,100*FullOutput2(2,11,2,2,1,1,2).Output.kappa);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). $\\kappa$ is the welfare benefit of introducing social security relative to a situation of no social security, measured by equivalent variation and expressed as a fraction of aggregate income. \n');
% fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). $\\kappa$ is the welfare benefit of introducing social security relative to a situation of no social security, measured by equivalent variation and expressed as a percentage of aggregate income. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% % Some outputs related to Table 3 that I used to check it was working
% [FullOutput2(2,2,2,2,1,1,2).Output.kappa,FullOutput2(2,3,2,2,1,1,2).Output.kappa,FullOutput2(2,4,2,2,1,1,2).Output.kappa,FullOutput2(2,5,2,2,1,1,2).Output.kappa,FullOutput2(2,6,2,2,1,1,2).Output.kappa,FullOutput2(2,7,2,2,1,1,2).Output.kappa,FullOutput2(2,11,2,2,1,1,2).Output.kappa]
% [FullOutput2(2,2,2,2,1,1,2).Output.minval,FullOutput2(2,3,2,2,1,1,2).Output.minval,FullOutput2(2,4,2,2,1,1,2).Output.minval,FullOutput2(2,5,2,2,1,1,2).Output.minval,FullOutput2(2,6,2,2,1,1,2).Output.minval,FullOutput2(2,7,2,2,1,1,2).Output.minval,FullOutput2(2,11,2,2,1,1,2).Output.minval]


% Table 4
FID = fopen('./SavedOutput/LatexInputs/ImrohorogluImrohorogluJoines1995_Table4.tex', 'w');
fprintf(FID, 'The role of intertemporal elasticity of substitution \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{2}{l}{$1/\\gamma=0.67$ ($\\gamma=1.5$)} & \\multicolumn{2}{l}{$1/\\gamma=0.50$ ($\\gamma=2.0$)} & \\multicolumn{2}{l}{$1/\\gamma=0.25$ ($\\gamma=4.0$)}  \\\\  \\cline{2-3} \\cline{4-5} \\cline{6-7} \n');
fprintf(FID, ' b & K/Q & Utility & K/Q & Utility & K/Q & Utility  \\\\ \\hline \n');
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,1,2,2,1,1,2).Output.Params.b,FullOutput2(2,1,2,1,1,1,2).Output.KdivQ,FullOutput2(2,1,2,1,1,1,2).Output.Omega1,FullOutput2(2,1,2,2,1,1,2).Output.KdivQ,FullOutput2(2,1,2,2,1,1,2).Output.Omega1,FullOutput2(2,1,2,3,1,1,2).Output.KdivQ,FullOutput2(2,1,2,3,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.Params.b,FullOutput2(2,2,2,1,1,1,2).Output.KdivQ,FullOutput2(2,2,2,1,1,1,2).Output.Omega1,FullOutput2(2,2,2,2,1,1,2).Output.KdivQ,FullOutput2(2,2,2,2,1,1,2).Output.Omega1,FullOutput2(2,2,2,3,1,1,2).Output.KdivQ,FullOutput2(2,2,2,3,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,3,2,2,1,1,2).Output.Params.b,FullOutput2(2,3,2,1,1,1,2).Output.KdivQ,FullOutput2(2,3,2,1,1,1,2).Output.Omega1,FullOutput2(2,3,2,2,1,1,2).Output.KdivQ,FullOutput2(2,3,2,2,1,1,2).Output.Omega1,FullOutput2(2,3,2,3,1,1,2).Output.KdivQ,FullOutput2(2,3,2,3,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,4,2,2,1,1,2).Output.Params.b,FullOutput2(2,4,2,1,1,1,2).Output.KdivQ,FullOutput2(2,4,2,1,1,1,2).Output.Omega1,FullOutput2(2,4,2,2,1,1,2).Output.KdivQ,FullOutput2(2,4,2,2,1,1,2).Output.Omega1,FullOutput2(2,4,2,3,1,1,2).Output.KdivQ,FullOutput2(2,4,2,3,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,5,2,2,1,1,2).Output.Params.b,FullOutput2(2,5,2,1,1,1,2).Output.KdivQ,FullOutput2(2,5,2,1,1,1,2).Output.Omega1,FullOutput2(2,5,2,2,1,1,2).Output.KdivQ,FullOutput2(2,5,2,2,1,1,2).Output.Omega1,FullOutput2(2,5,2,3,1,1,2).Output.KdivQ,FullOutput2(2,5,2,3,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,6,2,2,1,1,2).Output.Params.b,FullOutput2(2,6,2,1,1,1,2).Output.KdivQ,FullOutput2(2,6,2,1,1,1,2).Output.Omega1,FullOutput2(2,6,2,2,1,1,2).Output.KdivQ,FullOutput2(2,6,2,2,1,1,2).Output.Omega1,FullOutput2(2,6,2,3,1,1,2).Output.KdivQ,FullOutput2(2,6,2,3,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,7,2,2,1,1,2).Output.Params.b,FullOutput2(2,7,2,1,1,1,2).Output.KdivQ,FullOutput2(2,7,2,1,1,1,2).Output.Omega1,FullOutput2(2,7,2,2,1,1,2).Output.KdivQ,FullOutput2(2,7,2,2,1,1,2).Output.Omega1,FullOutput2(2,7,2,3,1,1,2).Output.KdivQ,FullOutput2(2,7,2,3,1,1,2).Output.Omega1);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Table 5
FID = fopen('./SavedOutput/LatexInputs/ImrohorogluImrohorogluJoines1995_Table5.tex', 'w');
fprintf(FID, 'The role of the discount factor and productivity growth \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llllllllll} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{3}{l}{$\\beta=0.98$,$g=0.0$} & \\multicolumn{3}{l}{$\\beta=1.011$,$g=0.022$} & \\multicolumn{3}{l}{$\\beta=1.011$,$g=0.0$}  \\\\  \\cline{2-4} \\cline{5-7} \\cline{8-10} \n');
fprintf(FID, ' b & K/Q & r & Utility & K/Q & r & Utility & K/Q & r & Utility  \\\\ \\hline \n');
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,1,2,2,1,1,2).Output.Params.b,FullOutput1(1,1,2,2,1,1,2).Output.KdivQ,FullOutput1(1,1,2,2,1,1,2).Output.r,FullOutput1(1,1,2,2,1,1,2).Output.Omega1,FullOutput2(2,1,2,2,2,1,2).Output.KdivQ,FullOutput2(2,1,2,2,2,1,2).Output.r,FullOutput2(2,1,2,2,2,1,2).Output.Omega1,FullOutput2(2,1,2,2,1,1,2).Output.KdivQ,FullOutput2(2,1,2,2,1,1,2).Output.r,FullOutput2(2,1,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.Params.b,FullOutput1(1,2,2,2,1,1,2).Output.KdivQ,FullOutput1(1,2,2,2,1,1,2).Output.r,FullOutput1(1,2,2,2,1,1,2).Output.Omega1,FullOutput2(2,2,2,2,2,1,2).Output.KdivQ,FullOutput2(2,2,2,2,2,1,2).Output.r,FullOutput2(2,2,2,2,2,1,2).Output.Omega1,FullOutput2(2,2,2,2,1,1,2).Output.KdivQ,FullOutput2(2,2,2,2,1,1,2).Output.r,FullOutput2(2,2,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,3,2,2,1,1,2).Output.Params.b,FullOutput1(1,3,2,2,1,1,2).Output.KdivQ,FullOutput1(1,3,2,2,1,1,2).Output.r,FullOutput1(1,3,2,2,1,1,2).Output.Omega1,FullOutput2(2,3,2,2,2,1,2).Output.KdivQ,FullOutput2(2,3,2,2,2,1,2).Output.r,FullOutput2(2,3,2,2,2,1,2).Output.Omega1,FullOutput2(2,3,2,2,1,1,2).Output.KdivQ,FullOutput2(2,3,2,2,1,1,2).Output.r,FullOutput2(2,3,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,4,2,2,1,1,2).Output.Params.b,FullOutput1(1,4,2,2,1,1,2).Output.KdivQ,FullOutput1(1,4,2,2,1,1,2).Output.r,FullOutput1(1,4,2,2,1,1,2).Output.Omega1,FullOutput2(2,4,2,2,2,1,2).Output.KdivQ,FullOutput2(2,4,2,2,2,1,2).Output.r,FullOutput2(2,4,2,2,2,1,2).Output.Omega1,FullOutput2(2,4,2,2,1,1,2).Output.KdivQ,FullOutput2(2,4,2,2,1,1,2).Output.r,FullOutput2(2,4,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,5,2,2,1,1,2).Output.Params.b,FullOutput1(1,5,2,2,1,1,2).Output.KdivQ,FullOutput1(1,5,2,2,1,1,2).Output.r,FullOutput1(1,5,2,2,1,1,2).Output.Omega1,FullOutput2(2,5,2,2,2,1,2).Output.KdivQ,FullOutput2(2,5,2,2,2,1,2).Output.r,FullOutput2(2,5,2,2,2,1,2).Output.Omega1,FullOutput2(2,5,2,2,1,1,2).Output.KdivQ,FullOutput2(2,5,2,2,1,1,2).Output.r,FullOutput2(2,5,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,6,2,2,1,1,2).Output.Params.b,FullOutput1(1,6,2,2,1,1,2).Output.KdivQ,FullOutput1(1,6,2,2,1,1,2).Output.r,FullOutput1(1,6,2,2,1,1,2).Output.Omega1,FullOutput2(2,6,2,2,2,1,2).Output.KdivQ,FullOutput2(2,6,2,2,2,1,2).Output.r,FullOutput2(2,6,2,2,2,1,2).Output.Omega1,FullOutput2(2,6,2,2,1,1,2).Output.KdivQ,FullOutput2(2,6,2,2,1,1,2).Output.r,FullOutput2(2,6,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,7,2,2,1,1,2).Output.Params.b,FullOutput1(1,7,2,2,1,1,2).Output.KdivQ,FullOutput1(1,7,2,2,1,1,2).Output.r,FullOutput1(1,7,2,2,1,1,2).Output.Omega1,FullOutput2(2,7,2,2,2,1,2).Output.KdivQ,FullOutput2(2,7,2,2,2,1,2).Output.r,FullOutput2(2,7,2,2,2,1,2).Output.Omega1,FullOutput2(2,7,2,2,1,1,2).Output.KdivQ,FullOutput2(2,7,2,2,1,1,2).Output.r,FullOutput2(2,7,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f & %8.3f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,11,2,2,1,1,2).Output.Params.b,FullOutput1(1,11,2,2,1,1,2).Output.KdivQ,FullOutput1(1,11,2,2,1,1,2).Output.r,FullOutput1(1,11,2,2,1,1,2).Output.Omega1,FullOutput2(2,11,2,2,2,1,2).Output.KdivQ,FullOutput2(2,11,2,2,2,1,2).Output.r,FullOutput2(2,11,2,2,2,1,2).Output.Omega1,FullOutput2(2,11,2,2,1,1,2).Output.KdivQ,FullOutput2(2,11,2,2,1,1,2).Output.r,FullOutput2(2,11,2,2,1,1,2).Output.Omega1);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


% Table 6
FID = fopen('./SavedOutput/LatexInputs/ImrohorogluImrohorogluJoines1995_Table6.tex', 'w');
fprintf(FID, 'Risk of catastrophic illness in old age \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllll} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{2}{l}{Prob. of illness 0.18} & \\multicolumn{2}{l}{Prob. of illness 0.09} & \\multicolumn{2}{l}{No catastrophic}  \\\\ \n');
fprintf(FID, ' & \\multicolumn{2}{l}{cost 25\\%%} & \\multicolumn{2}{l}{cost 35\\%%} & \\multicolumn{2}{l}{illness}  \\\\  \\cline{2-3} \\cline{4-5} \\cline{6-7} \n');
fprintf(FID, ' b & K/Q & Utility & K/Q & Utility & K/Q & Utility  \\\\ \\hline \n');
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,1,2,2,1,1,2).Output.Params.b,FullOutput2(2,1,2,2,1,2,2).Output.KdivQ,FullOutput2(2,1,2,2,1,2,2).Output.Omega1,FullOutput2(2,1,2,2,1,3,2).Output.KdivQ,FullOutput2(2,1,2,2,1,3,2).Output.Omega1,FullOutput2(2,1,2,2,1,1,2).Output.KdivQ,FullOutput2(2,1,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,2,2,2,1,1,2).Output.Params.b,FullOutput2(2,2,2,2,1,2,2).Output.KdivQ,FullOutput2(2,2,2,2,1,2,2).Output.Omega1,FullOutput2(2,2,2,2,1,3,2).Output.KdivQ,FullOutput2(2,2,2,2,1,3,2).Output.Omega1,FullOutput2(2,2,2,2,1,1,2).Output.KdivQ,FullOutput2(2,2,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,3,2,2,1,1,2).Output.Params.b,FullOutput2(2,3,2,2,1,2,2).Output.KdivQ,FullOutput2(2,3,2,2,1,2,2).Output.Omega1,FullOutput2(2,3,2,2,1,3,2).Output.KdivQ,FullOutput2(2,3,2,2,1,3,2).Output.Omega1,FullOutput2(2,3,2,2,1,1,2).Output.KdivQ,FullOutput2(2,3,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,4,2,2,1,1,2).Output.Params.b,FullOutput2(2,4,2,2,1,2,2).Output.KdivQ,FullOutput2(2,4,2,2,1,2,2).Output.Omega1,FullOutput2(2,4,2,2,1,3,2).Output.KdivQ,FullOutput2(2,4,2,2,1,3,2).Output.Omega1,FullOutput2(2,4,2,2,1,1,2).Output.KdivQ,FullOutput2(2,4,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,5,2,2,1,1,2).Output.Params.b,FullOutput2(2,5,2,2,1,2,2).Output.KdivQ,FullOutput2(2,5,2,2,1,2,2).Output.Omega1,FullOutput2(2,5,2,2,1,3,2).Output.KdivQ,FullOutput2(2,5,2,2,1,3,2).Output.Omega1,FullOutput2(2,5,2,2,1,1,2).Output.KdivQ,FullOutput2(2,5,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,6,2,2,1,1,2).Output.Params.b,FullOutput2(2,6,2,2,1,2,2).Output.KdivQ,FullOutput2(2,6,2,2,1,2,2).Output.Omega1,FullOutput2(2,6,2,2,1,3,2).Output.KdivQ,FullOutput2(2,6,2,2,1,3,2).Output.Omega1,FullOutput2(2,6,2,2,1,1,2).Output.KdivQ,FullOutput2(2,6,2,2,1,1,2).Output.Omega1);
fprintf(FID, ' %8.2f & %8.3f & %8.2f & %8.3f & %8.2f & %8.3f & %8.2f  \\\\ \n',FullOutput2(2,7,2,2,1,1,2).Output.Params.b,FullOutput2(2,7,2,2,1,2,2).Output.KdivQ,FullOutput2(2,7,2,2,1,2,2).Output.Omega1,FullOutput2(2,7,2,2,1,3,2).Output.KdivQ,FullOutput2(2,7,2,2,1,3,2).Output.Omega1,FullOutput2(2,7,2,2,1,1,2).Output.KdivQ,FullOutput2(2,7,2,2,1,1,2).Output.Omega1);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: b is the social security replacement rate (IIJ1995 call this $\\theta$). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


