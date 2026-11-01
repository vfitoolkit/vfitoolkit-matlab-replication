function Out=AS2000_SolveColumn(Params,Grids,vfoptions,simoptions,heteroagentoptions,doGE)
% Solves one column of Tables 2-5: the stationary equilibrium for given parameters.
% doGE=1: find the tax rate tau that balances the UI budget (start from Params.tau,
%         a low guess, as the paper looks for the lowest balanced-budget tau, fn 12)
% doGE=0: tau is taken as given (No UI column: uib=0 and tau=0)

n_d=Grids.n_d;
d_grid=Grids.d_grid;
n_a=Grids.n_a;
a_grid=Grids.a_grid;
n_semiz=vfoptions.n_semiz;
% No exogenous z, all the risk is in the semi-exogenous employment state
n_z=0;
z_grid=[];
pi_z=[];

DiscountFactorParamNames={'beta','Gammadisc'}; % beta*Gamma^(1-theta), eq (5)

ReturnFn=@(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,R,Gamma,hbar,eta,theta,ojs) AS2000_ReturnFn(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,R,Gamma,hbar,eta,theta,ojs);

%% General eqm: government budget balance
if doGE==1
    FnsToEvaluateGE.TaxRevenue=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) tau*AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2));
    FnsToEvaluateGE.UIspending=@(h,s,aprime,a,semiz,uib) uib*(semiz==3);
    GeneralEqmEqns.govbudget=@(TaxRevenue,UIspending) TaxRevenue-UIspending;
    GEPriceParamNames={'tau'};

    p_eqm=HeteroAgentStationaryEqm_InfHorz(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluateGE,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    Params.tau=p_eqm.tau;
end

%% Solve at the equilibrium tau
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params,[]);

%% Statistics
% Hours: work hours when employed, search hours when unemployed
FnsToEvaluate.assets=@(h,s,aprime,a,semiz) a;
FnsToEvaluate.good=@(h,s,aprime,a,semiz) (semiz==1);
FnsToEvaluate.bad=@(h,s,aprime,a,semiz) (semiz==2);
FnsToEvaluate.insured=@(h,s,aprime,a,semiz) (semiz==3);
FnsToEvaluate.uninsured=@(h,s,aprime,a,semiz) (semiz==4);
FnsToEvaluate.hoursgood=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(semiz==1);
FnsToEvaluate.hoursbad=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(semiz==2);
FnsToEvaluate.hoursinsured=@(h,s,aprime,a,semiz) h*(semiz==3);
FnsToEvaluate.hoursuninsured=@(h,s,aprime,a,semiz) h*(semiz==4);
FnsToEvaluate.searchbad=@(h,s,aprime,a,semiz,ojs) ojs*h*(semiz==2); % on-the-job search hours (Table 6)
% Weekly (pre-tax) labor income. The std dev of its log is over the employed with positive earnings (with on-the-job
% search a bad-job worker can be all-search, with zero earnings); the log uses log(1)=0 for the rest, who are multiplied by zero
FnsToEvaluate.earnings=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2));
FnsToEvaluate.posearnings=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) (AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))>0);
FnsToEvaluate.logearnings=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) (AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))>0)*log(AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))+(AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))<=0));
FnsToEvaluate.logearnings2=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) (AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))>0)*log(AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))+(AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))<=0))^2;
FnsToEvaluate.wage=@(h,s,aprime,a,semiz,wg,wb) wg*(semiz==1)+wb*(semiz==2);
% Leisure is reported as 1-h, not hbar-h: Table 3 (hbar=1.7) reports leisure 0.150 with hours ~0.86,
% which only fits 1-h (with hbar=1, as in Tables 2,4,5, the two are the same)
FnsToEvaluate.leisure=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) 1-AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)-ojs*h*(semiz==2);
% Flow into unemployment: quits (s=1 when employed), plus exogenous separations of those who stay
FnsToEvaluate.inflowU=@(h,s,aprime,a,semiz,deltag,deltab) (semiz==1)*((s==1)+(s==0)*deltag)+(semiz==2)*((s==1)+(s==0)*deltab);
FnsToEvaluate.quit=@(h,s,aprime,a,semiz) (semiz<=2)*(s==1);
FnsToEvaluate.seekgood=@(h,s,aprime,a,semiz) (semiz>=3)*(s==1);
FnsToEvaluate.consumption=@(h,s,aprime,a,semiz,wg,wb,uib,uincfloor,tau,R,Gamma,hbar,eta,ojs) a+(1-tau)*AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))+max(uib,uincfloor)*(semiz==3)+uincfloor*(semiz==4)-aprime*Gamma/R;
FnsToEvaluate.TaxRevenue=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) tau*AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2));
FnsToEvaluate.UIspending=@(h,s,aprime,a,semiz,uib) uib*(semiz==3);
FnsToEvaluate.hours=@(h,s,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs) AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs);
FnsToEvaluate.s=@(h,s,aprime,a,semiz) s;

AggVars=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

E=AggVars.good.Mean+AggVars.bad.Mean; % employment
U=AggVars.insured.Mean+AggVars.uninsured.Mean; % unemployment

Out.uib=Params.uib;
Out.phi=Params.phi;
Out.tau=Params.tau;
Out.meanassets=AggVars.assets.Mean;
Out.fracgood=AggVars.good.Mean/E;
Out.urate=U;
Out.uduration=U/AggVars.inflowU.Mean; % steady state: mean completed duration = stock/inflow
Eposearn=AggVars.posearnings.Mean; % employed with positive earnings
Out.sdlogweeklyincome=sqrt(AggVars.logearnings2.Mean/Eposearn-(AggVars.logearnings.Mean/Eposearn)^2);
Out.fraczeroearnings=1-Eposearn/E; % share of the employed with zero earnings (all-search, only with on-the-job search)
Out.hoursgood=AggVars.hoursgood.Mean/AggVars.good.Mean;
Out.hoursbad=AggVars.hoursbad.Mean/AggVars.bad.Mean;
Out.searchbad=AggVars.searchbad.Mean/AggVars.bad.Mean;
Out.hoursinsured=AggVars.hoursinsured.Mean/AggVars.insured.Mean;
Out.hoursuninsured=AggVars.hoursuninsured.Mean/AggVars.uninsured.Mean;
Out.hoursallunemployed=(AggVars.hoursinsured.Mean+AggVars.hoursuninsured.Mean)/U; % for the No UI column (insured=uninsured when uib=0)
Out.avgwage=AggVars.wage.Mean/E;
Out.output=AggVars.earnings.Mean;
Out.leisure=AggVars.leisure.Mean;
Out.consumption=AggVars.consumption.Mean;
Out.quitrate=AggVars.quit.Mean/E;
Out.fracseekgood=AggVars.seekgood.Mean/U;
Out.govbudget=AggVars.TaxRevenue.Mean-AggVars.UIspending.Mean;
Out.UIspending=AggVars.UIspending.Mean;

% Welfare, consumption metric: the c* that, consumed forever without working (leisure hbar),
% gives the same utility as the average worker in the stationary dist.
betatilde=Params.beta*Params.Gammadisc;
V=gather(V);
StationaryDist=gather(StationaryDist);
EV=sum(V(StationaryDist>0).*StationaryDist(StationaryDist>0));
Out.EV=EV;
if Params.theta==1
    Out.welfare=exp((1-betatilde)*EV-Params.eta*log(Params.hbar));
else
    Out.welfare=((1-Params.theta)*(1-betatilde)*EV)^(1/(1-Params.theta))/(Params.hbar^Params.eta);
end

% Mass at the top of the asset grid, and search hours at the edges of the h grid (grid checks)
Out.massattopofagrid=sum(StationaryDist(end,:));
Out.massabove90pctagrid=sum(sum(StationaryDist(a_grid>0.9*a_grid(end),:)));
hgrid=d_grid(1:n_d(1));
hoursU=reshape(gather(ValuesOnGrid.hours),[n_a,prod(n_semiz)]);
hoursU=hoursU(:,3:4);
distU=StationaryDist(:,3:4);
Out.massUathgridbottom=sum(distU(hoursU<=hgrid(1)+1e-12))/U;
Out.massUathgridtop=sum(distU(hoursU>=hgrid(end)-1e-12))/U;

%% Simulated panel: N=1000 workers for T=10000 weeks, starting from the stationary dist
% Used for the consumption drop following job loss, the std dev of log annual income (both pooled over
% the T and N dimensions), and Fig 6 (the asset path of the first worker).
FnsToEvaluatePanel.assets=@(h,s,aprime,a,semiz) a;
FnsToEvaluatePanel.consumption=FnsToEvaluate.consumption;
% Income Y_t (Section 3): after-tax labor income when employed, benefits when insured, zero when uninsured
% (the 1e-7 income floor is not counted as income)
FnsToEvaluatePanel.income=@(h,s,aprime,a,semiz,wg,wb,uib,tau,R,Gamma,hbar,eta,ojs) (1-tau)*AS2000_Hours(h,aprime,a,semiz,wg,wb,tau,R,Gamma,hbar,eta,ojs)*(wg*(semiz==1)+wb*(semiz==2))+uib*(semiz==3);
FnsToEvaluatePanel.semiz=@(h,s,aprime,a,semiz) semiz;
simoptionsPanel=simoptions;
simoptionsPanel.numbersims=1000;
simoptionsPanel.simperiods=10000;
simoptionsPanel.burnin=0; % InitialDist is the stationary dist
SimPanel=SimPanelValues_InfHorz(StationaryDist,Policy,FnsToEvaluatePanel,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptionsPanel); % each field is [T,N]
panelc=gather(SimPanel.consumption);
panely=gather(SimPanel.income);
panelz=gather(SimPanel.semiz);

% Consumption drop following job loss: every week t where an employed worker (good or bad job) is laid off,
% ie. becomes insured unemployed in t+1 (quits, which make the worker uninsured, are not job loss).
% Drop is 1-c(t+1)/c(t), averaged over all these events.
joblossevent=(panelz(1:end-1,:)<=2) & (panelz(2:end,:)==3);
consdrop=1-panelc(2:end,:)./panelc(1:end-1,:);
Out.consdrop=mean(consdrop(joblossevent));
Out.njoblossevents=sum(joblossevent(:));
% Std dev of log annual income: sum income over 52-week years (192 complete years per worker),
% pooled over all worker-years. Years with zero income (uninsured all year) have no log, and are counted.
nyears=floor(simoptionsPanel.simperiods/52);
annualincome=reshape(sum(reshape(panely(1:52*nyears,:),[52,nyears,simoptionsPanel.numbersims]),1),[nyears*simoptionsPanel.numbersims,1]);
Out.sdlogannualincome=std(log(annualincome(annualincome>0)),1);
Out.nzeroincomeyears=sum(annualincome==0);
% Fig 6: asset path of the first worker
Out.assetpath=gather(SimPanel.assets(:,1));
% Cross-check of the panel against the stationary dist
Out.panelmeanassets=mean(mean(gather(SimPanel.assets)));

% Policies on grid, [n_a,n_semiz]
Out.cpolicy=reshape(gather(ValuesOnGrid.consumption),[n_a,prod(n_semiz)]);
Out.hpolicy=reshape(gather(ValuesOnGrid.hours),[n_a,prod(n_semiz)]);
Out.spolicy=reshape(gather(ValuesOnGrid.s),[n_a,prod(n_semiz)]);
Out.StationaryDist=reshape(StationaryDist,[n_a,prod(n_semiz)]);
Out.V=reshape(V,[n_a,prod(n_semiz)]);
Out.Policy=gather(Policy);
Out.Params=Params;

end
