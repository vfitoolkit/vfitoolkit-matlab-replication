% RunBaseline_gridconv.m  --  DIAGNOSTIC (not the deliverable).
%
% The var(log c) RISE (our headline replication statistic) moved ~0.025 when amax went
% 60x->100x at a fixed 381-point positive segment, revealing it is NOT grid-converged at
% n_a=501 (W/Y and delta* ARE converged; the rise is not; see NOTES.md). Two effects were
% confounded there: curing top-truncation (lowers the rise) and coarsening the populated
% region (unclear). This isolates BULK resolution: fix amax=100x (so no truncation) and
% delta=0.9495 (the calibrated lambda=0.62 value; W/Y is flat near it so it stays ~4), and
% sweep ONLY the number of positive-segment asset points until the rise stabilises.
% lambda=0.62 only. If the rise is flat by the top of the sweep, that value is the
% converged headline number; if still moving, the deliverable grid must be refined further.

clearvars -except doPart
if exist('RunBaseline_gridconv_diary.txt','file'); delete('RunBaseline_gridconv_diary.txt'); end
diary('RunBaseline_gridconv_diary.txt');
fprintf('=== RunBaseline_gridconv.m (lambda=0.62, amax=100x, delta=0.9495 fixed) run %s ===\n',char(datetime('now')));

%% Parameters (identical to RunBaseline; lambda=0.62, delta fixed at its calibrated value)
Jwork=40; Jret=31; N_j=Jwork+Jret;
s2a_u=0.022; s2b_u=0.00038; sab_u=-0.00045;
Params.rho=0.821; Params.sigma2_eta=0.029; Params.sigma2_eps=0.047;
Params.betabar=0.009; Params.alphabar=1.5; Params.exper1=0; Params.T=Jwork;
R2ab=sab_u^2/(s2a_u*s2b_u);
Params.sigma2_alpha=0; Params.sigma_ab=0; Params.sigma2_beta=s2b_u*(1-R2ab);
Params.crra=2; Params.Pb=0.96; Params.r=1/Params.Pb-1;
Params.delta=0.9495;         % calibrated lambda=0.62 value (W/Y flat near it -> stays ~4)
Params.Phibar=1/1.40;
Params.exper_j=[0:Jwork-1,(Jwork-1)*ones(1,Jret)];
Params.exper_T=Jwork-1;
Params.workret_j=[ones(1,Jwork),zeros(1,Jret)];
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};
DiscountFactorParamNames={'delta'};
lambda=0.62; Params.lambda=lambda;

n_d=0; d_grid=[];
n_z=[61,13,5]; N_z=prod(n_z); nb=n_z(1);
Params.nSigmasBelief=3;
discopts.parallel=0; discopts.nSigmas=2.5;
n_types=5;

% positive-segment point counts to sweep (negative segment fixed at 121)
npos_list=[61, 121, 181]; % totals 181/241/301 (capped for GI memory at N_z=3965)

%% alpha types
[alphagrid,w_types]=discretizeIIDNormal_TanakaToda(Params.alphabar,sqrt(s2a_u),n_types,discopts);
alphagrid=gather(alphagrid); w_types=gather(w_types);
dbet=(sab_u/s2a_u)*(alphagrid-Params.alphabar);

%% Natural borrowing limit
ymin=zeros(1,N_j); zmin=0;
for jj=1:Jwork
    zmin=Params.rho*zmin-2.5*sqrt(Params.sigma2_eta);
    ymin(jj)=exp((Params.alphabar-2.5*sqrt(s2a_u))+(Params.betabar-2.5*sqrt(s2b_u))*Params.exper_j(jj)+zmin-2.5*sqrt(Params.sigma2_eps));
end
ymin(Jwork+1:N_j)=Params.Phibar*0.9*ymin(Jwork);
M=zeros(1,N_j+1); for jj=N_j:-1:1; M(jj)=ymin(jj)+Params.Pb*M(jj+1); end
Params.Wbar_j=[-M(2:N_j),0];

%% Belief chain (built once; only the wealth grid changes)
kf=KalmanSetup(Params);
[z_grid_J,pi_z_J,jequaloneDistz,grids]=BuildPiZ(kf,Params,n_z,N_j,discopts);
pj=jequaloneDistz'; for jj=1:Jwork-1; pj=pj*pi_z_J(:,:,jj); end
ynode_T=grids.bgrid_J(grids.IB,Jwork)*Params.exper_T+grids.zgrid_J(grids.IZ,Jwork)+grids.vgrid_J(grids.IV,Jwork);
Params.Ybar_T=(pj*exp(ynode_T))*(w_types'*exp(alphagrid+dbet*Params.exper_T));
amin=min(Params.Wbar_j); amax=100*Params.Ybar_T;
fprintf('Ybar_T=%.3f, amax=%.1f, delta=%.4f\n',Params.Ybar_T,amax,Params.delta);

%% FnsToEvaluate
FnsToEvaluate.cons=@(aprime,a,bh,zh,v,workret_j,exper_j,Pb,Phibar,Ybar_T,exper_T,alpha_i) ...
    a-Pb*aprime+workret_j*exp(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<0.3)*(0.9*exp(alpha_i+bh*exper_T+zh+v)/Ybar_T) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>=0.3)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=2)*(0.27+0.32*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-0.3)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>2)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-2)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>4.1)*1.1);
FnsToEvaluate.income=@(aprime,a,bh,zh,v,workret_j,exper_j,Phibar,Ybar_T,exper_T,alpha_i) ...
    workret_j*exp(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<0.3)*(0.9*exp(alpha_i+bh*exper_T+zh+v)/Ybar_T) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>=0.3)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=2)*(0.27+0.32*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-0.3)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>2)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-2)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>4.1)*1.1);
FnsToEvaluate.assets=@(aprime,a,bh,zh,v) a;
FnsToEvaluate.logc=@(aprime,a,bh,zh,v,workret_j,exper_j,Pb,Phibar,Ybar_T,exper_T,alpha_i) ...
    log(max(1e-12,a-Pb*aprime+workret_j*exp(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<0.3)*(0.9*exp(alpha_i+bh*exper_T+zh+v)/Ybar_T) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>=0.3)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=2)*(0.27+0.32*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-0.3)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>2)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-2)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>4.1)*1.1)));
ReturnFn=@(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i) ...
    ReturnFn_HIP(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i);
vfoptions.verbose=0; vfoptions.divideandconquer=1;
vfoptions.gridinterplayer=1;  % interpolate aprime between asset nodes (smooth policy)
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=1; % MUST match vfoptions for StationaryDist/LifeCycleProfiles
simoptions.ngridinterp=10;


%% Sweep positive-segment resolution
ns=numel(npos_list);
res_na=zeros(1,ns); res_rise=zeros(1,ns); res_v65=zeros(1,ns); res_WY=zeros(1,ns); res_top=zeros(1,ns); res_t=zeros(1,ns);
for is=1:ns
    npos=npos_list(is);
    a_grid=unique([linspace(amin,0,121)'; (amax*(linspace(0,1,npos).^3))']);
    n_a=length(a_grid);
    [~,a0ind]=min(abs(a_grid));
    jequaloneDist=zeros(n_a,N_z,'gpuArray'); jequaloneDist(a0ind,:)=jequaloneDistz';
    jequaloneDist=reshape(jequaloneDist,[n_a,n_z]);
    Ma=zeros(n_types,N_j); My=zeros(n_types,N_j); Mlc=zeros(n_types,N_j); Vlc=zeros(n_types,N_j);
    topmass=0; tic;
    for ii=1:n_types
        Params.alpha_i=alphagrid(ii);
        zgJ=z_grid_J; zgJ(1:nb,:)=z_grid_J(1:nb,:)+dbet(ii);
        [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        SD=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
        ACS=LifeCycleProfiles_FHorz_Case1(SD,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ,simoptions);
        Ma(ii,:)=gather(ACS.assets.Mean); My(ii,:)=gather(ACS.income.Mean);
        Mlc(ii,:)=gather(ACS.logc.Mean); Vlc(ii,:)=gather(ACS.logc.Variance);
        topmass=topmass+w_types(ii)*gather(sum(SD(end-2:end,:,:,:,:),'all'));
    end
    vft=toc;
    meana=w_types'*Ma; meany=w_types'*My;
    varlogc=w_types'*(Vlc+Mlc.^2)-(w_types'*Mlc).^2;
    res_na(is)=n_a; res_rise(is)=varlogc(Jwork+1)-varlogc(1); res_v65(is)=varlogc(Jwork+1);
    res_WY(is)=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
    res_top(is)=topmass; res_t(is)=vft;
    fprintf('npos=%4d (n_a=%4d): rise=%.4f, var(logc)_65=%.4f, W/Y=%.2f, top-mass=%.1e, %.0fs\n',...
        npos,n_a,res_rise(is),res_v65(is),res_WY(is),topmass,vft);
end

fprintf('\n--- Convergence of the var(log c) rise (lambda=0.62) ---\n');
fprintf('n_a:      ');fprintf(' %8d',res_na);fprintf('\n');
fprintf('rise:     ');fprintf(' %8.4f',res_rise);fprintf('\n');
fprintf('W/Y:      ');fprintf(' %8.2f',res_WY);fprintf('\n');
fprintf('successive rise changes: ');fprintf(' %8.4f',diff(res_rise));fprintf('\n');
fprintf('READING: if the last successive change is < ~0.005 the rise is converged; take the\n');
fprintf('  finest npos as the deliverable positive-segment size. If still moving, refine further.\n');

save('RunBaseline_gridconv_results.mat','npos_list','res_na','res_rise','res_v65','res_WY','res_top','res_t','Params');
fprintf('\n=== RunBaseline_gridconv done ===\n');
diary off;
