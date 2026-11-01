% Acemoglu & Shimer (2000) - Productivity gains from unemployment insurance
% European Economic Review 44, 1195-1224.
%
% Replicates the dynamic model (Sections 3-5): Tables 2-6, Figs 3-6 and the Section 4.5 numbers. Rows 13-14 of the
% tables and Fig 6 use a simulated panel of N=1000 workers for T=10000 weeks (see AS2000_SolveColumn.m).
% (Figs 1-2 are the static model of Section 2, not replicated.)
%
% Infinite horizon, one period is a week. Workers save in a riskless bond (a>=0) at exogenous R.
% Everything is normalized by productivity Gamma^t: a'=(R/Gamma)*(a+y-c), discount factor beta*Gamma^(1-theta).
% Employment state is a semi-exogenous state (its transitions depend on decisions):
%   semiz: 1=good job, 2=bad job, 3=unemployed insured, 4=unemployed uninsured
% Decisions (both are d2, ie. both affect the semiz transitions, vfoptions.l_dsemiz=2):
%   h: search hours (only matters when unemployed; employed hours come from the static FOC inside the ReturnFn)
%   s: unemployed: s=1 seek a good job, s=0 seek a bad job; employed: s=1 quit, s=0 stay
% Transitions (decided this period, take effect next period):
%   employed who stay lose the job with prob delta -> insured unemployed
%   employed who quit -> uninsured unemployed (we allow quitting good jobs too; no-one should)
%   unemployed find a job of the type they seek with prob m*h; insured who do not find a job lose eligibility with prob phi
% General eqm: tax rate tau on labor income balances the UI budget, tau*E[h*w]=b*mass(insured).
% The paper takes the lowest balanced-budget tau (fn 12), so we start from a low guess.

%% Run infrastructure (diary, paths, output folders)
diary off
if exist('./AS2000_Diary.txt','file'), delete('./AS2000_Diary.txt'); end
diary ./AS2000_Diary.txt

close all

addpath('./AS2000subcodes');
% To run on server, we need to tell Matlab where to find VFI Toolkit
if exist('./VFIToolkit-matlab','dir'); addpath(genpath('./VFIToolkit-matlab/')); end

if ~exist('./SavedOutput','dir'), mkdir('./SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'), mkdir('./SavedOutput/Graphs'); end

doPart=[1,1,1,1,1,1];
% doPart(1): Table 2 (benchmark parameterization, four UI schemes), Figs 3-6, and the Section 4.5 numbers
% doPart(2): Table 3 (high value of leisure, eta=1, hbar=1.7)
% doPart(3): Table 4 (low wage dispersion, wg=1.15, mb=0.15, mg=0.06)
% doPart(4): Table 5 (low risk aversion, theta=1, wg=1.15, mb=0.2, mg=0.04)
% doPart(5): Table 6 (on-the-job search in bad jobs)
% doPart(6): the LaTeX tables of the write-up (./SavedOutput/LatexInputs/), from the saved results; solves nothing
% Each doPart saves its results to ./SavedOutput/AS2000_doPartX.mat

%% Parameters, Table 1
% Technology
Params.R=1.0006;      % gross weekly interest rate
Params.Gamma=1.0003;  % gross weekly productivity growth rate
Params.wb=1;          % bad job wage
Params.deltab=0.005;  % bad job layoff rate
Params.mb=0.5;        % bad job matching rate
Params.wg=1.3;        % good job wage
Params.deltag=0.005;  % good job layoff rate
Params.mg=0.05;       % good job matching rate
% Preferences
Params.beta=0.998;    % discount factor
Params.theta=4;       % relative risk aversion
Params.hbar=1;        % maximum hours
Params.eta=1/6;       % relative value of leisure (Table 1 reports 0.167, text says 1/6)
% UI policy
Params.uib=0.4;       % unemployment benefit
Params.phi=0.04;      % benefit exhaustion rate
% Tiny income of the unemployed (not in the paper). Without it an unemployed worker with zero assets has
% no feasible consumption (c=0 gives u=-Inf), so V(a=0,uninsured)=-Inf (and V(a=0,insured)=-Inf in the
% No UI column). In the paper's continuous-asset model this state is never reached from a>0, but on a grid
% it puts a -Inf into the value function, and the expected value next to it can turn into NaN. With this
% floor the state is merely very bad (u is about -7e20 per period), which keeps everything finite and,
% being 1e-7 of a bad-job wage, has no effect on the results. The insured get max(uib,uincfloor), so the floor
% only bites for them when uib=0 (No UI). It is not counted in output or UI spending.
Params.uincfloor=10^(-7);
% On-the-job search (Section 5.4): =1 lets bad-job workers search for good jobs, used for Table 6 only
Params.ojs=0;
% General eqm
Params.tau=0.022;     % tax rate (initial guess)

Params.Gammadisc=Params.Gamma^(1-Params.theta); % growth-adjustment to the discount factor, eq (5)

%% Grids
% Assets: the paper uses {1.01^n-1} for about 400 points (fn 15); no-one holds more than 37 in the benchmark
n_a=396;
a_grid=1.01.^(0:n_a-1)'-1; % max is 50.6
% Decisions: h (search hours) and s (binary)
n_h=101;
h_grid=linspace(0.5,1,n_h)'; % search hours, observed range in Tables 2,4,5 is 0.78-0.92
s_grid=[0;1];
n_d=[n_h,2];
d_grid=[h_grid; s_grid];

Grids.n_a=n_a;
Grids.a_grid=a_grid;
Grids.n_d=n_d;
Grids.d_grid=d_grid;

%% Toolkit options
vfoptions.n_semiz=4;
vfoptions.semiz_grid=(1:1:4)';
vfoptions.l_dsemiz=2; % both decisions affect the semiz transitions
vfoptions.SemiExoStateFn=@(semiz,semizprime,h,s,deltag,deltab,mg,mb,phi,ojs) AS2000_SemiExoStateFn(semiz,semizprime,h,s,deltag,deltab,mg,mb,phi,ojs);
vfoptions.gridinterplayer=1; % a' is chosen from a finer grid (paper interpolates the value fn linearly); without it a' gets stuck at a
vfoptions.ngridinterp=20;
vfoptions.tolerance=1e-9;
vfoptions.verbose=0;

simoptions.n_semiz=vfoptions.n_semiz;
simoptions.semiz_grid=vfoptions.semiz_grid;
simoptions.l_dsemiz=vfoptions.l_dsemiz;
simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions.d_grid=d_grid;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;
simoptions.tolerance=1e-10; % asset holdings are very persistent (Fig 6), so the dist converges slowly
simoptions.verbose=0;

heteroagentoptions.verbose=1;

%% The paper's numbers (rows as in Table 2, columns Benchmark, High UI, Long UI, No UI; NaN for a dash)
PaperTable2=[0.4,0.44,0.4,0; 25,25,27.8,NaN; 0.022,0.025,0.023,0; 16.5,16.8,17.0,7.1;...
    51.5,56.1,54.2,0; 6.3,6.7,6.5,1.1; 13.5,14.5,14.0,2.2; 0.133,0.132,0.132,0.002;...
    0.862,0.861,0.862,NaN; 0.859,0.859,0.859,0.858; 0.816,0.809,0.814,NaN; 0.866,0.866,0.866,0.895;...
    15.2,15.4,15.4,10.2; 0.169,0.168,0.171,0.020;...
    1.155,1.169,1.163,1.000; 0.931,0.938,0.935,0.849; 0.142,0.142,0.142,0.142; 0.651,0.655,0.654,0.611];
PaperTable3=[0.4,0.44,0.4,0; 25,25,27.8,NaN; 0.025,0.031,0.027,0; 16.2,15.5,16.1,7.1;...
    53.0,55.4,53.8,0; 7.3,8.0,7.7,1.0; 15.8,17.5,16.6,2.0; 0.138,0.138,0.138,0.008;...
    0.866,0.865,0.866,NaN; 0.855,0.854,0.855,0.853; 0.672,0.624,0.663,NaN; 0.886,0.893,0.888,0.971;...
    22.4,24.2,24.0,4.2; 0.181,0.183,0.183,0.014;...
    1.160,1.168,1.162,1; 0.926,0.923,0.924,0.844; 0.150,0.155,0.151,0.146; 0.446,0.446,0.446,0.419];
PaperTable4=[0.4,0.44,0.4,0; 25,25,27.8,NaN; 0.023,0.026,0.024,0; 27.1,26.9,27.1,27.1;...
    46.5,49.2,48.8,0; 6.3,6.5,6.4,3.7; 13.5,14.0,13.7,7.6; 0.070,0.070,0.070,0.005;...
    0.860,0.860,0.860,NaN; 0.859,0.859,0.859,0.860; 0.822,0.816,0.820,NaN; 0.875,0.878,0.875,0.873;...
    15.0,15.3,15.1,12.2; 0.136,0.136,0.136,0.133;...
    1.070,1.074,1.073,1; 0.862,0.863,0.863,0.828; 0.142,0.143,0.142,0.140; 0.608,0.609,0.610,0.588];
PaperTable5=[0.4,0.44,0.4,0; 25,25,27.8,NaN; 0.025,0.033,0.028,0; 5.8,6.0,5.9,5.8;...
    33.8,46.7,38.3,0; 6.4,7.7,6.9,2.6; 13.7,16.8,14.7,5.4; 0.066,0.070,0.068,0.004;...
    0.859,0.859,0.859,NaN; 0.859,0.859,0.859,0.859; 0.828,0.821,0.827,NaN; 0.921,0.921,0.922,0.910;...
    27.7,27.9,27.9,28.6; 0.114,0.134,0.121,0.068;...
    1.051,1.070,1.057,1; 0.845,0.848,0.846,0.836; 0.142,0.143,0.142,0.140; 0.603,0.604,0.604,0.598];
% Table 6 has the extra row 'Search hours in bad jobs' after 'Work hours in bad jobs' (19 rows)
% (Leisure in the No UI column is printed as 0.862 in the paper, presumably a typo)
PaperTable6=[0.4,0.44,0.4,0; 25,25,27.8,NaN; 0.021,0.025,0.022,0; 17.8,17.5,18.0,15.0;...
    61.1,63.1,62.9,37.4; 6.4,6.7,6.6,3.2; 13.7,14.5,14.2,6.6; 0.212,0.208,0.168,0.199;...
    0.862,0.861,0.862,0.863; 0.839,0.840,0.839,0.835; 0.020,0.021,0.022,0.027;...
    0.818,0.810,0.816,NaN; 0.865,0.865,0.865,0.864;...
    15.6,15.9,16.7,13.0; 0.172,0.173,0.173,0.176;...
    1.185,1.191,1.190,1.115; 0.946,0.948,0.949,0.912; 0.142,0.141,0.138,0.862; 0.659,0.660,0.660,0.637];

%% Tables 2-6
TableNames={'Table 2: benchmark parameterization','Table 3: high value of leisure, eta=1',...
    'Table 4: low wage dispersion, wg=1.15','Table 5: low risk aversion, theta=1','Table 6: on-the-job search'};
PaperTables={PaperTable2,PaperTable3,PaperTable4,PaperTable5,PaperTable6};
ColNames={'Benchmark','High UI','Long UI','No UI'};

for tt=1:5
    if doPart(tt)==1
        fprintf('doPart(%i): %s \n',tt,TableNames{tt})
        tic1=tic;
        ParamsT=Params;
        GridsT=Grids;
        simoptionsT=simoptions;
        if tt==2 % Section 5.1
            ParamsT.eta=1;
            ParamsT.hbar=1.7;
            % Search hours of the insured fall to about 0.62 (Table 3), and hbar is higher
            h_gridT=linspace(0.4,1.4,n_h)';
            GridsT.n_d=[length(h_gridT),2];
            GridsT.d_grid=[h_gridT; s_grid];
            simoptionsT.d_grid=GridsT.d_grid;
        elseif tt==3 % Section 5.2
            ParamsT.wg=1.15;
            ParamsT.mb=0.15;
            ParamsT.mg=0.06;
        elseif tt==4 % Section 5.3
            ParamsT.theta=1;
            ParamsT.wg=1.15;
            ParamsT.mb=0.2;
            ParamsT.mg=0.04;
        elseif tt==5 % Section 5.4
            % Bad-job workers choose work hours hw and search hours hs (the h decision), finding a good job with
            % prob mg*hs (if the bad job does not end first). Leisure is hbar-hw-hs. All other parameters as Table 2.
            ParamsT.ojs=1;
            % Search hours on the job are mostly zero, and all-search (about the same as an unemployed worker)
            % when assets pass a threshold (Section 5.4), so the h grid now starts at zero
            h_gridT=[linspace(0,0.5,51)'; h_grid(2:end)];
            GridsT.n_d=[length(h_gridT),2];
            GridsT.d_grid=[h_gridT; s_grid];
            simoptionsT.d_grid=GridsT.d_grid;
        end
        ParamsT.Gammadisc=ParamsT.Gamma^(1-ParamsT.theta);

        Outs=cell(1,4);
        for cc=1:4
            ParamsC=ParamsT;
            ParamsC.tau=PaperTables{tt}(3,cc); % initial guess: the paper's (lowest) tau
            doGE=1;
            if cc==2 % High UI: benefits up 10%
                ParamsC.uib=0.44;
            elseif cc==3 % Long UI: exhaustion rate down 10%
                ParamsC.phi=0.9*ParamsT.phi;
            elseif cc==4 % No UI
                ParamsC.uib=0;
                ParamsC.tau=0;
                doGE=0;
            end
            tic2=tic;
            Outs{cc}=AS2000_SolveColumn(ParamsC,GridsT,vfoptions,simoptionsT,heteroagentoptions,doGE);
            O=Outs{cc};
            fprintf('  %s: tau=%.5f (budget residual %.2e), took %.1f minutes \n',ColNames{cc},O.tau,O.govbudget,toc(tic2)/60)
            fprintf('    grid checks: mass at top of a_grid %.2e, mass above 0.9*max(a_grid) %.2e, unemployed at bottom/top of h_grid %.2e / %.2e \n',...
                O.massattopofagrid,O.massabove90pctagrid,O.massUathgridbottom,O.massUathgridtop)
            % Thresholds (paper, benchmark: insured seek good jobs above about 4, uninsured above about 9, bad jobs quit above 19)
            for zz=[3,4]
                switches=sum(abs(diff(O.spolicy(:,zz)))>0);
                ind=find(O.spolicy(:,zz)==1,1,'first');
                if isempty(ind)
                    fprintf('    semiz=%i: never seeks good jobs \n',zz)
                else
                    fprintf('    semiz=%i: seeks good jobs from a=%.2f (policy switches %i times) \n',zz,a_grid(ind),switches)
                end
            end
            for zz=[1,2]
                ind=find(O.spolicy(:,zz)==1,1,'first');
                if isempty(ind)
                    fprintf('    semiz=%i: never quits \n',zz)
                else
                    fprintf('    semiz=%i: quits from a=%.2f \n',zz,a_grid(ind))
                end
            end
            fprintf('    quit rate (share of employed) %.2e, share of unemployed seeking good jobs %.3f \n',O.quitrate,O.fracseekgood)
            fprintf('    panel: %i job-loss events, %i zero-income years (dropped from log annual income), mean assets %.3f (stationary dist %.3f) \n',...
                O.njoblossevents,O.nzeroincomeyears,O.panelmeanassets,O.meanassets)
            if ParamsC.ojs==1
                fprintf('    on-the-job search: mean search hours in bad jobs %.3f, share of employed with zero earnings (all-search) %.4f \n',O.searchbad,O.fraczeroearnings)
            end
        end
        AS2000_PrintTable(TableNames{tt},Outs,PaperTables{tt},1);

        %% Section 4.5 (Table 2 only): size of the gains, and the bound on the transition costs
        if tt==1
            O1=Outs{1};
            fprintf('Section 4.5: gains from UI relative to the benchmark (paper in brackets) \n')
            PaperS45=[17.6, 0.0035, 0.7, 0.6; 5.6, 0.0011, 0.5, 0.4]; % UI spending % change, UI spending change per worker, output % change, welfare % change
            for cc=2:3
                O=Outs{cc};
                dUI=O.UIspending-O1.UIspending;
                fprintf('  %s: UI spending %+.1f%% (%.1f%%), %+.4f per worker (%.4f); output %+.2f%% (%.1f%%); welfare %+.2f%% (%.1f%%); change in output / change in UI spending %.2f \n',...
                    ColNames{cc},100*dUI/O1.UIspending,PaperS45(cc-1,1),dUI,PaperS45(cc-1,2),100*(O.output/O1.output-1),PaperS45(cc-1,3),100*(O.welfare/O1.welfare-1),PaperS45(cc-1,4),(O.output-O1.output)/dUI)
            end
            % Transition costs: start from the benchmark stationary dist, and permanently switch to the High UI (Long UI)
            % policy together with its tax rate. The worker's problem does not depend on the distribution, so their value
            % fn is exactly that of the High UI (Long UI) column, and welfare is that V integrated over the benchmark dist.
            % (The government runs a surplus early on, which is not redistributed, so this is an upper bound on the transition costs.)
            fprintf('Section 4.5: welfare including transition, from the benchmark stationary dist (paper: 0.651 to 0.652 for both) \n')
            S45.welfaretransition=zeros(1,3);
            for cc=2:3
                O=Outs{cc};
                Pc=O.Params;
                betatilde=Pc.beta*Pc.Gammadisc;
                EVtrans=sum(O.V(O1.StationaryDist>0).*O1.StationaryDist(O1.StationaryDist>0));
                if Pc.theta==1
                    S45.welfaretransition(cc)=exp((1-betatilde)*EVtrans-Pc.eta*log(Pc.hbar));
                else
                    S45.welfaretransition(cc)=((1-Pc.theta)*(1-betatilde)*EVtrans)^(1/(1-Pc.theta))/(Pc.hbar^Pc.eta);
                end
                fprintf('  %s: benchmark %.4f, with transition %.4f, ignoring transition (Table 2) %.4f \n',ColNames{cc},O1.welfare,S45.welfaretransition(cc),O.welfare)
            end
            save('./SavedOutput/AS2000_Section45.mat','S45')
        end
        fprintf('doPart(%i) took %.1f minutes \n',tt,toc(tic1)/60)
        save(['./SavedOutput/AS2000_doPart',num2str(tt),'.mat'],'Outs','ParamsT','GridsT')

        %% Figs 3-6 (benchmark column of Table 2)
        if tt==1
            O=Outs{1};
            a_grid=GridsT.a_grid;
            plotind=(a_grid<=40);
            % Bad-job workers who quit are not in a bad job, so do not plot their policy there
            cplot=O.cpolicy; hplot=O.hpolicy;
            cplot(O.spolicy(:,2)==1,2)=NaN; hplot(O.spolicy(:,2)==1,2)=NaN;
            cplot(O.spolicy(:,1)==1,1)=NaN; hplot(O.spolicy(:,1)==1,1)=NaN;
            PanelOrder=[4,3,2,1];
            PanelTitles={'Unemployed & Uninsured','Unemployed & Insured','Employed in Bad Job','Employed in Good Job'};

            fig3=figure(3);
            for pp=1:4
                subplot(2,2,pp); plot(a_grid(plotind),cplot(plotind,PanelOrder(pp)));
                title(PanelTitles{pp}); xlabel('normalized assets');
            end
            sgtitle('Fig 3: consumption')
            saveas(fig3,'./SavedOutput/Graphs/AS2000_Fig3.png')

            fig4=figure(4);
            for pp=1:4
                subplot(2,2,pp); plot(a_grid(plotind),hplot(plotind,PanelOrder(pp)));
                title(PanelTitles{pp}); xlabel('normalized assets');
            end
            sgtitle('Fig 4: hours of search or labor supply')
            saveas(fig4,'./SavedOutput/Graphs/AS2000_Fig4.png')

            % Density of assets: mass at each grid point divided by the width of the grid cell it represents
            massa=sum(O.StationaryDist,2);
            da=gradient(a_grid);
            fig5=figure(5);
            plot(a_grid(plotind),massa(plotind)./da(plotind));
            title('Fig 5: ergodic distribution of asset holdings'); xlabel('normalized assets');
            saveas(fig5,'./SavedOutput/Graphs/AS2000_Fig5.png')

            % Same distribution as a cdf (the pdf on a grid is sensitive to numerical error, the cdf is not)
            fig51=figure(51);
            plot(a_grid(plotind),cumsum(massa(plotind)));
            title('Fig 5 (cdf): ergodic distribution of asset holdings'); xlabel('normalized assets');
            saveas(fig51,'./SavedOutput/Graphs/AS2000_Fig5cdf.png')

            % Asset path of the first worker in the simulated panel
            fig6=figure(6);
            plot(1:length(O.assetpath),O.assetpath);
            title('Fig 6: a 10,000 week sample path of asset holdings'); xlabel('week');
            saveas(fig6,'./SavedOutput/Graphs/AS2000_Fig6.png')
        end
    end
end

%% doPart(6): LaTeX tables for the write-up
if doPart(6)==1
    AS2000_LatexTables('./SavedOutput');
end

diary off
