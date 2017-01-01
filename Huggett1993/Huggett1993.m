% Replicates the results (and more) of Huggett (1993)
% Figures use plotly
% Note: The algorithm used for calculating the general equilibrium is not what you would want to use normally for solving this model. It finds the general equilibrium using a discrete grid on interest rates, rather
% than just solving the fixed-point problem on interest rates directly by using optimization. This is done for robustness reasons and because these codes were written as part of a paper on
% a convergent algorithm for solving Bewley-Huggett-Aiyagari-models.


Parallel=2 %GPU


% A few lines needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
try % Server has 16 cores, but is shared with other users, so use max of 8.
    parpool(8)
catch % Desktop has less than 8, so will give error, on desktop it is fine to use all available cores.
    parpool
end
PoolDetails=gcp;
NCores=PoolDetails.NumWorkers;

%% Set some basic variables

n_a=2^10;% Huggett used "between 150 and 350 evenly spaced gridpoints" (and uses linear interpolation)
n_e=2;
n_q=551;

%Parameters
Params.beta=0.99322; %Model period is one-sixth of a year, so implied 'annual beta' is 0.96.

Params.mu=1.5; % Huggett (1993) calls this sigma
mu_vec=[1.5,3];

Params.alowerbar=-2; % {-2,-4,-6,-8}
alowerbar_vec=[-2,-4,-6,-8];
%Params.q % The price

Params.eh=1;
Params.el=0.1;
Params.pieheh=0.925; % Probability of eh given eh
Params.piehel=0.5; % Probability of eh given el

%% Some Toolkit options
simoptions.ncores=NCores; % Number of CPU cores
vfoptions=struct(); % Just use default options

%% Solve the model a bunch of times
Tables=struct();
TimeTable=zeros(4,2);
for alowerbar_c=1:4
    Params.alowerbar=alowerbar_vec(alowerbar_c);
    for mu_c=1:2
        Params.mu=mu_vec(mu_c);
        fprintf('Current iteration alowerbar_c=%d, mu_c=%d \n', alowerbar_c,mu_c)
        tic;
        Output=Huggett1993_Fn(n_a,n_e,n_q,Params, Parallel, vfoptions, simoptions);
        time1=toc
        % Output is a structure
        Tables(alowerbar_c,mu_c).q=Output.q;
        Tables(alowerbar_c,mu_c).r=Output.r;
        
        Figures(alowerbar_c,mu_c).Policy=Output.Policy;
        Figures(alowerbar_c,mu_c).StationaryDist=Output.StationaryDist;
        
        Figures(alowerbar_c,mu_c).a_grid=Output.a_grid;
        
        TimeTable(alowerbar_c,mu_c)=time1;
    end
end

save ./SavedOutput/Huggett1993Tables.mat Tables TimeTable

%% Reproduce Tables 1 & 2 of Huggett (1993)

mu_c=1;
FilenameString=['./SavedOutput/LatexInputs/Huggett1993_Table1.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, 'Coefficient of Relative Risk Aversion $\\mu$=%1.1f \\\\ \n', mu_vec(mu_c));
fprintf(FID, '\\begin{tabular*}{0.5\\textwidth}{@{\\extracolsep{\\fill}}lcr} \\hline \\hline \n');
fprintf(FID, ' Credit Limit & Interest Rate & Price \\\\ \n');
fprintf(FID, ' (-$\\underbar{a}$) & ($r$) & ($q$) \\\\ \\hline \n');
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(1), Tables(1,mu_c).r, Tables(1,mu_c).q);
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(2), Tables(2,mu_c).r, Tables(2,mu_c).q);
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(3), Tables(3,mu_c).r, Tables(3,mu_c).q);
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(4), Tables(4,mu_c).r, Tables(4,mu_c).q);
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Replication of Table 1 of Huggett (1993) using grid sizes $ n_a=%d $, $ n_e=%d $, $ n_q=%d $ \n', n_a, n_e, n_q);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

mu_c=2;
FilenameString=['./SavedOutput/LatexInputs/Huggett1993_Table2.tex'];
FID = fopen(FilenameString, 'w');
fprintf(FID, 'Coefficient of Relative Risk Aversion $\\mu$=%1.1f \\\\ \n', mu_vec(mu_c));
fprintf(FID, '\\begin{tabular*}{0.5\\textwidth}{@{\\extracolsep{\\fill}}lcr} \\hline \\hline \n');
fprintf(FID, ' Credit Limit & Interest Rate & Price \\\\ \n');
fprintf(FID, ' (-$\\underbar{a}$) & ($r$) & ($q$) \\\\ \\hline \n');
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(1), Tables(1,mu_c).r, Tables(1,mu_c).q);
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(2), Tables(2,mu_c).r, Tables(2,mu_c).q);
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(3), Tables(3,mu_c).r, Tables(3,mu_c).q);
fprintf(FID, '  %d & %1.1f \\%% & %1.4f \\\\ \n', alowerbar_vec(4), Tables(4,mu_c).r, Tables(4,mu_c).q);
fprintf(FID, '\\hline \\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Replication of Table 2 of Huggett (1993) using grid sizes $ n_a=%d $, $ n_e=%d $, $ n_q=%d $ \n', n_a, n_e, n_q);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Reproduce Figures 1 & 2 of Huggett (1993)
% % Note that figures depend on plotly
% a_grid=Figures(1,1).a_grid;
% 
% % Figure 1
% Policy=Figures(1,1).Policy;
% trace1= struct('x', a_grid,'y', Policy(:,1),'name','e_l','type', 'scatter');
% trace2= struct('x', a_grid,'y', Policy(:,2),'name','e_h','type', 'scatter');
% trace3= struct('x', a_grid,'y', a_grid,'name','45 deg line','type', 'scatter');
% data = {trace1,trace2, trace3};
% layout = struct('title', 'Optimal Decision Rule','showlegend', true,'width', 800,...
%     'xaxis', struct('domain', [0, 1],'title','Assets (a)','showgrid',false), ...
%     'yaxis', struct('title', 'Decision (next period assets)','showgrid',false,'titlefont', struct('color', 'black'),'tickfont', struct('color', 'black'),'side', 'left','position',0),...
%     'legend',struct('x',0.8,'y',0.8));
% response = plotly(data, struct('layout', layout, 'filename', 'Huggett1993_Figure1', 'fileopt', 'overwrite'));
% response.data=data; response.layout=layout;
% saveplotlyfig(response, ['./SavedOutput/Graphs/Huggett1993_Figure1.pdf'])
% 
% % Figure 2
% StationaryDist=Figures(1,1).StationaryDist;
% trace1= struct('x', a_grid,'y', cumsum(StationaryDist(:,1)),'name','e_l','type', 'scatter');
% trace2= struct('x', a_grid,'y', cumsum(StationaryDist(:,2)),'name','e_h','type', 'scatter');
% data = {trace1,trace2};
% layout = struct('title', 'Stationarty Distribution','showlegend', true,'width', 800,...
%     'xaxis', struct('domain', [0, 1],'title','Assets (a)','showgrid',false), ...
%     'yaxis', struct('title', 'Cumulative Distribution Fn','showgrid',false,'titlefont', struct('color', 'black'),'tickfont', struct('color', 'black'),'side', 'left','position',0),...
%     'legend',struct('x',0.8,'y',0.8));
% response = plotly(data, struct('layout', layout, 'filename', 'Huggett1993_Figure2', 'fileopt', 'overwrite'));
% response.data=data; response.layout=layout;
% saveplotlyfig(response, ['./SavedOutput/Huggett1993_Figure2.pdf'])



