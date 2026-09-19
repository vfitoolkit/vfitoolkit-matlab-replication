% RunEducationGroups.m
% Stage (d), Section III.C / Figure 11: the co-movement of consumption and income by
% EDUCATION group. Solves the college and high-school calibrations (Table 1 rows 3-6)
% for both HIP and RIP, each delta-calibrated to W/Y=4 on the accurate grid, and reports
% the average life-cycle consumption GROWTH by group. The paper's discriminating result:
%   HIP: consumption growth DIFFERS by education (college +53% vs HS +29%, ages 25-55)
%   RIP: consumption growth is NEARLY IDENTICAL (43.4% vs 42.7%) -- counterfactual.
% This is an aggregate mean-consumption object; our distribution method gives it directly
% (no panel simulation). HIP uses the known-alpha reduction (5 alpha types, n_z=[61,13,5]);
% RIP uses alpha as an identity-transition z-dimension (n_z=[7,15,5]), as in RunRIP.
% Requires GPU. Run from baseline/.

clearvars -except doPart
if exist('RunEducationGroups_diary.txt','file'); delete('RunEducationGroups_diary.txt'); end
if ~exist('SavedOutput/Graphs','dir'); mkdir('SavedOutput/Graphs'); end
diary('RunEducationGroups_diary.txt');
fprintf('=== RunEducationGroups.m run %s ===\n',char(datetime('now')));

%% ============================ SAB CHECK (up front) ============================
% Table 1 prints sigma_ab per group, but for the all-sample we established the printed
% value is rescaled: the results require sab=-0.00045, not the printed -0.0020 (factor
% 0.225; reproduces Table 2, Figs 3/4). There is NO per-education Table 2 to cross-check
% the groups, so we cannot pin them the same way. We report both readings and default to
% the same rescale, because it yields COHERENT modest correlations across all rows
% (all -0.156, college -0.161, HS -0.057) consistent with the gradual-learning story,
% whereas the printed values imply implausibly strong correlations (-0.69/-0.72/-0.25).
% NB: sab is SECOND-ORDER for the Figure-11 consumption-GROWTH statistic (driven mainly by
% meanbeta, s2b and lambda); it mainly shifts the between-type tilt. Flag, don't agonise.
sab_rule='rescaled';   % 'rescaled' (default) or 'printed'
resc=0.00045/0.0020;   % all-sample true/printed = 0.225
% [s2a, s2b, sab_printed] per HIP group (Table 1 rows 4,6)
grp_s2a =[0.023, 0.038];   % college, HS
grp_s2b =[0.00049, 0.00020];
grp_sabP=[-0.0024, -0.0007];
fprintf('\n--- SAB CHECK (rule=%s) ---\n',sab_rule);
fprintf('%-10s %12s %12s %12s %12s\n','group','sab_printed','corr_print','sab_used','corr_used');
grp_sab=zeros(1,2);
nm={'college','highschl'};
for g=1:2
    if strcmp(sab_rule,'printed')
        grp_sab(g)=grp_sabP(g); %#ok<UNRCH> % live only when the flag above is set to 'printed'
    else
        grp_sab(g)=grp_sabP(g)*resc;
    end
    cP=grp_sabP(g)/sqrt(grp_s2a(g)*grp_s2b(g));
    cU=grp_sab(g)/sqrt(grp_s2a(g)*grp_s2b(g));
    fprintf('%-10s %12.5f %12.3f %12.6f %12.3f\n',nm{g},grp_sabP(g),cP,grp_sab(g),cU);
end
fprintf('(second-order for the consumption-growth statistic; set sab_rule to test sensitivity)\n');

%% ============================ common calibration ============================
Jwork=40; Jret=31; N_j=Jwork+Jret; agevec=25:95;
Params.crra=2; Params.Pb=0.96; Params.r=1/Params.Pb-1;
Params.Phibar=1/1.40; Params.alphabar=1.5; Params.exper1=0; Params.T=Jwork;
Params.exper_j=[0:Jwork-1,(Jwork-1)*ones(1,Jret)];
Params.exper_T=Jwork-1;
Params.workret_j=[ones(1,Jwork),zeros(1,Jret)];
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};
DiscountFactorParamNames={'delta'};
n_d=0; d_grid=[];
discopts.parallel=0; discopts.nSigmas=2.5;
vfoptions.verbose=0;
vfoptions.gridinterplayer=1;  % interpolate aprime between asset nodes (smooth policy)
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=1; % MUST match vfoptions for StationaryDist/LifeCycleProfiles
simoptions.ngridinterp=10;

% group income-process parameters (Table 1 rows 3-6; meanbeta from Table 3)
% HIP: rows 4 (college), 6 (HS); RIP: rows 3 (college), 5 (HS)
edu(1)=struct('name','college','model','HIP');
edu(2)=struct('name','highschl','model','HIP');
edu(3)=struct('name','college','model','RIP');
edu(4)=struct('name','highschl','model','RIP');
% HIP params
hipP(1)=struct('rho',0.805,'s2a',0.023,'s2b',0.00049,'sab',grp_sab(1),'s2eta',0.025,'s2eps',0.032,'lambda',0.55,'betabar',0.012); % college
hipP(2)=struct('rho',0.829,'s2a',0.038,'s2b',0.00020,'sab',grp_sab(2),'s2eta',0.022,'s2eps',0.034,'lambda',0.32,'betabar',0.007); % HS
% RIP params (rows 3,5); betabar per education
ripP(1)=struct('rho',0.979,'s2a',0.031,'s2eta',0.0099,'s2eps',0.047,'betabar',0.012); % college
ripP(2)=struct('rho',0.972,'s2a',0.053,'s2eta',0.011,'s2eps',0.052,'betabar',0.007); % HS

results=struct();

%% ============================ HIP groups ============================
n_z=[61,13,5]; N_z=prod(n_z); nb=n_z(1); Params.nSigmasBelief=3;
n_types=5;
for g=1:2
    P=hipP(g);
    fprintf('\n============ HIP %s (lambda=%.2f, betabar=%.3f) ============\n',edu(g).name,P.lambda,P.betabar);
    R2ab=P.sab^2/(P.s2a*P.s2b);
    Params.rho=P.rho; Params.sigma2_eta=P.s2eta; Params.sigma2_eps=P.s2eps;
    Params.betabar=P.betabar; Params.lambda=P.lambda;
    Params.sigma2_alpha=0; Params.sigma_ab=0; Params.sigma2_beta=P.s2b*(1-R2ab); % conditional-on-alpha
    s2a_u=P.s2a; s2b_u=P.s2b; sab_u=P.sab;

    % alpha types
    [alphagrid,w_types]=discretizeIIDNormal_TanakaToda(Params.alphabar,sqrt(s2a_u),n_types,discopts);
    alphagrid=gather(alphagrid); w_types=gather(w_types);
    dbet=(sab_u/s2a_u)*(alphagrid-Params.alphabar);

    % natural borrowing limit (this group's public-prior HIP params, unconditional s2b)
    ymin=zeros(1,N_j); zmin=0;
    for jj=1:Jwork
        zmin=Params.rho*zmin-2.5*sqrt(Params.sigma2_eta);
        ymin(jj)=exp((Params.alphabar-2.5*sqrt(s2a_u))+(Params.betabar-2.5*sqrt(s2b_u))*Params.exper_j(jj)+zmin-2.5*sqrt(Params.sigma2_eps));
    end
    ymin(Jwork+1:N_j)=Params.Phibar*0.9*ymin(Jwork);
    M=zeros(1,N_j+1); for jj=N_j:-1:1; M(jj)=ymin(jj)+Params.Pb*M(jj+1); end
    Params.Wbar_j=[-M(2:N_j),0];

    % belief chain + Ybar_T + asset grid
    kf=KalmanSetup(Params);
    [z_grid_J,pi_z_J,jequaloneDistz,grids]=BuildPiZ(kf,Params,n_z,N_j,discopts);
    pj=jequaloneDistz'; for jj=1:Jwork-1; pj=pj*pi_z_J(:,:,jj); end
    ynode_T=grids.bgrid_J(grids.IB,Jwork)*Params.exper_T+grids.zgrid_J(grids.IZ,Jwork)+grids.vgrid_J(grids.IV,Jwork);
    Params.Ybar_T=(pj*exp(ynode_T))*(w_types'*exp(alphagrid+dbet*Params.exper_T));
    amin=min(Params.Wbar_j); amax=100*Params.Ybar_T;
    a_grid=unique([linspace(amin,0,97)'; (amax*(linspace(0,1,279).^3))']);
    n_a=length(a_grid);
    fprintf('R2=%.4f, cond s2b=%.6f, Ybar_T=%.3f, n_a=%d on [%.2f,%.2f]\n',R2ab,Params.sigma2_beta,Params.Ybar_T,n_a,amin,amax);

    % FnsToEvaluate (alpha_i per-type parameter; log guard)
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
    FnsCalib.assets=FnsToEvaluate.assets; FnsCalib.income=FnsToEvaluate.income;
    ReturnFn=@(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i) ...
        ReturnFn_HIP(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i);
    vfoptions.divideandconquer=1;

    [~,a0ind]=min(abs(a_grid));
    jequaloneDist=zeros(n_a,N_z,'gpuArray'); jequaloneDist(a0ind,:)=jequaloneDistz';
    jequaloneDist=reshape(jequaloneDist,[n_a,n_z]);

    % delta calibration to W/Y=4 (aggregate over types)
    dlo=0.90; dhi=0.978;
    for calibit=1:8
        Params.delta=(dlo+dhi)/2;
        meana_agg=zeros(1,N_j); meany_agg=zeros(1,N_j);
        for ii=1:n_types
            Params.alpha_i=alphagrid(ii);
            zgJ=z_grid_J; zgJ(1:nb,:)=z_grid_J(1:nb,:)+dbet(ii);
            [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
            SD=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
            ACS=LifeCycleProfiles_FHorz_Case1(SD,Policy,FnsCalib,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ,simoptions);
            meana_agg=meana_agg+w_types(ii)*gather(ACS.assets.Mean);
            meany_agg=meany_agg+w_types(ii)*gather(ACS.income.Mean);
        end
        WYtot=sum(Params.mewj.*meana_agg)/sum(Params.mewj.*(meany_agg+Params.r*meana_agg));
        if WYtot>4; dhi=Params.delta; else; dlo=Params.delta; end
    end
    Params.delta=(dlo+dhi)/2;

    % final solve + aggregate
    Mc=zeros(n_types,N_j); My=zeros(n_types,N_j); Ma=zeros(n_types,N_j);
    Mlc=zeros(n_types,N_j); Vlc=zeros(n_types,N_j);
    for ii=1:n_types
        Params.alpha_i=alphagrid(ii);
        zgJ=z_grid_J; zgJ(1:nb,:)=z_grid_J(1:nb,:)+dbet(ii);
        [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        SD=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
        ACS=LifeCycleProfiles_FHorz_Case1(SD,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ,simoptions);
        Mc(ii,:)=gather(ACS.cons.Mean); My(ii,:)=gather(ACS.income.Mean); Ma(ii,:)=gather(ACS.assets.Mean);
        Mlc(ii,:)=gather(ACS.logc.Mean); Vlc(ii,:)=gather(ACS.logc.Variance);
    end
    meanc=w_types'*Mc; meany=w_types'*My; meana=w_types'*Ma;
    varlogc=w_types'*(Vlc+Mlc.^2)-(w_types'*Mlc).^2;
    WYtot=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
    results(g).name=edu(g).name; results(g).model='HIP'; results(g).delta=Params.delta;
    results(g).meanc=meanc; results(g).meany=meany; results(g).varlogc=varlogc; results(g).WYtot=WYtot;
    results(g).cg55=meanc(31)/meanc(1)-1; results(g).cg64=meanc(40)/meanc(1)-1;
    fprintf('delta*=%.4f, W/Y=%.2f, cons growth 25-55=%.1f%% (25-64=%.1f%%), var(logc) rise=%.3f\n',...
        Params.delta,WYtot,100*results(g).cg55,100*results(g).cg64,varlogc(Jwork+1)-varlogc(1));
end

%% ============================ RIP groups ============================
n_alpha=7; n_zshock=15; n_eps=5; n_zR=[n_alpha,n_zshock,n_eps]; N_zR=prod(n_zR);
for g=3:4
    P=ripP(g-2);
    fprintf('\n============ RIP %s (betabar=%.3f) ============\n',edu(g).name,P.betabar);
    Params.rho=P.rho; Params.sigma2_alpha=P.s2a; Params.sigma2_eta=P.s2eta; Params.sigma2_eps=P.s2eps;
    Params.betabar=P.betabar;
    % borrowing limit uses this group's HIP public-prior params (comparability, p.702)
    H=hipP(g-2);
    [alpha_grid,w_alpha]=discretizeIIDNormal_TanakaToda(Params.alphabar,sqrt(Params.sigma2_alpha),n_alpha,discopts);
    [eps_grid,w_eps]=discretizeIIDNormal_TanakaToda(0,sqrt(Params.sigma2_eps),n_eps,discopts);
    [zshock_grid_J,pi_zshock_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(...
        zeros(1,Jwork),Params.rho*ones(1,Jwork),sqrt(Params.sigma2_eta)*ones(1,Jwork),n_zshock,Jwork,discopts);
    alpha_grid=gather(alpha_grid); w_alpha=gather(w_alpha); eps_grid=gather(eps_grid); w_eps=gather(w_eps);
    zshock_grid_J=gather(zshock_grid_J); pi_zshock_J=gather(pi_zshock_J); jequaloneDistz=gather(jequaloneDistz);
    z_grid_J=zeros(sum(n_zR),N_j);
    for jj=1:N_j; jz=min(jj,Jwork); z_grid_J(:,jj)=[alpha_grid; zshock_grid_J(:,jz); eps_grid]; end
    pi_eps_iid=ones(n_eps,1)*w_eps';
    pi_z_J=zeros(N_zR,N_zR,N_j-1);
    for jj=1:N_j-1
        if jj<=Jwork-1; pi_z_J(:,:,jj)=kron(pi_eps_iid,kron(pi_zshock_J(:,:,jj),eye(n_alpha)));
        else; pi_z_J(:,:,jj)=eye(N_zR); end
    end
    jointgrid_T=[repmat(alpha_grid,n_zshock*n_eps,1), ...
        repmat(kron(zshock_grid_J(:,Jwork),ones(n_alpha,1)),n_eps,1), kron(eps_grid,ones(n_alpha*n_zshock,1))];
    pjoint=kron(w_eps,kron(jequaloneDistz,w_alpha))';
    for jj=1:Jwork-1; pjoint=pjoint*pi_z_J(:,:,jj); end
    Params.Ybar_T=pjoint*exp(jointgrid_T(:,1)+Params.betabar*Params.exper_T+jointgrid_T(:,2)+jointgrid_T(:,3));
    % natural borrowing limit (group HIP public-prior params)
    ymin=zeros(1,N_j); zmin=0;
    for jj=1:Jwork
        zmin=H.rho*zmin-2.5*sqrt(H.s2eta);
        ymin(jj)=exp((Params.alphabar-2.5*sqrt(H.s2a))+(Params.betabar-2.5*sqrt(H.s2b))*Params.exper_j(jj)+zmin-2.5*sqrt(H.s2eps));
    end
    ymin(Jwork+1:N_j)=Params.Phibar*0.9*ymin(Jwork);
    M=zeros(1,N_j+1); for jj=N_j:-1:1; M(jj)=ymin(jj)+Params.Pb*M(jj+1); end
    Params.Wbar_j=[-M(2:N_j),0];
    amin=min(Params.Wbar_j); amax=100*Params.Ybar_T;
    a_grid=unique([linspace(amin,0,97)'; (amax*(linspace(0,1,279).^3))']); n_a=length(a_grid);
    fprintf('Ybar_T=%.3f, n_a=%d on [%.2f,%.2f]\n',Params.Ybar_T,n_a,amin,amax);

    FnsToEvaluate=struct();
    FnsToEvaluate.cons=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,Pb,Phibar,Ybar_T,exper_T) ...
        a-Pb*aprime+workret_j*exp(alpha+betabar*exper_j+z+e)+(1-workret_j)*Phibar*Ybar_T*( ...
        (exp(alpha+betabar*exper_T+z+e)/Ybar_T<0.3)*(0.9*exp(alpha+betabar*exper_T+z+e)/Ybar_T) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>=0.3)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=2)*(0.27+0.32*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-0.3)) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>2)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-2)) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>4.1)*1.1);
    FnsToEvaluate.income=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,Phibar,Ybar_T,exper_T) ...
        workret_j*exp(alpha+betabar*exper_j+z+e)+(1-workret_j)*Phibar*Ybar_T*( ...
        (exp(alpha+betabar*exper_T+z+e)/Ybar_T<0.3)*(0.9*exp(alpha+betabar*exper_T+z+e)/Ybar_T) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>=0.3)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=2)*(0.27+0.32*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-0.3)) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>2)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-2)) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>4.1)*1.1);
    FnsToEvaluate.assets=@(aprime,a,alpha,z,e) a;
    FnsToEvaluate.logc=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,Pb,Phibar,Ybar_T,exper_T) ...
        log(max(1e-12,a-Pb*aprime+workret_j*exp(alpha+betabar*exper_j+z+e)+(1-workret_j)*Phibar*Ybar_T*( ...
        (exp(alpha+betabar*exper_T+z+e)/Ybar_T<0.3)*(0.9*exp(alpha+betabar*exper_T+z+e)/Ybar_T) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>=0.3)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=2)*(0.27+0.32*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-0.3)) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>2)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-2)) ...
        +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>4.1)*1.1)));
    FnsCalib=struct(); FnsCalib.assets=FnsToEvaluate.assets; FnsCalib.income=FnsToEvaluate.income;
    ReturnFn=@(aprime,a,alpha,z,e,workret_j,exper_j,Wbar_j,betabar,Pb,crra,Phibar,Ybar_T,exper_T) ...
        ReturnFn_RIP(aprime,a,alpha,z,e,workret_j,exper_j,Wbar_j,betabar,Pb,crra,Phibar,Ybar_T,exper_T);
    vfoptions.divideandconquer=1;

    [~,a0ind]=min(abs(a_grid));
    jequaloneDist=zeros([n_a,n_zR],'gpuArray');
    jequaloneDist(a0ind,:,:,:)=reshape(kron(w_eps,kron(jequaloneDistz,w_alpha)),[1,n_zR]);

    dlo=0.90; dhi=0.978;
    for calibit=1:8
        Params.delta=(dlo+dhi)/2;
        [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zR,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        SD=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_zR,N_j,pi_z_J,Params,simoptions);
        ACS=LifeCycleProfiles_FHorz_Case1(SD,Policy,FnsCalib,Params,[],n_d,n_a,n_zR,N_j,d_grid,a_grid,z_grid_J,simoptions);
        meana=gather(ACS.assets.Mean); meany=gather(ACS.income.Mean);
        WYtot=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
        if WYtot>4; dhi=Params.delta; else; dlo=Params.delta; end
    end
    Params.delta=(dlo+dhi)/2;
    [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zR,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    SD=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_zR,N_j,pi_z_J,Params,simoptions);
    ACS=LifeCycleProfiles_FHorz_Case1(SD,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_zR,N_j,d_grid,a_grid,z_grid_J,simoptions);
    meanc=gather(ACS.cons.Mean); meany=gather(ACS.income.Mean); meana=gather(ACS.assets.Mean);
    varlogc=gather(ACS.logc.Variance);
    WYtot=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
    results(g).name=edu(g).name; results(g).model='RIP'; results(g).delta=Params.delta;
    results(g).meanc=meanc; results(g).meany=meany; results(g).varlogc=varlogc; results(g).WYtot=WYtot;
    results(g).cg55=meanc(31)/meanc(1)-1; results(g).cg64=meanc(40)/meanc(1)-1;
    fprintf('delta*=%.4f, W/Y=%.2f, cons growth 25-55=%.1f%% (25-64=%.1f%%), var(logc) rise=%.3f\n',...
        Params.delta,WYtot,100*results(g).cg55,100*results(g).cg64,varlogc(Jwork+1)-varlogc(1));
end

%% ============================ summary: Figure 11 test ============================
fprintf('\n=== Figure 11: consumption growth by education (ages 25-55), HIP vs RIP ===\n');
fprintf('%-8s %-8s %10s %10s %10s\n','model','group','cons25-55','cons25-64','delta*');
for g=1:4
    fprintf('%-8s %-8s %9.1f%% %9.1f%% %10.4f\n',results(g).model,results(g).name,...
        100*results(g).cg55,100*results(g).cg64,results(g).delta);
end
fprintf('\nDISCRIMINATING TEST:\n');
fprintf('  HIP college - HS growth gap (25-55): %.1f pp (paper: 53%%-29%%=24pp -> DIFFERS)\n',...
    100*(results(1).cg55-results(2).cg55));
fprintf('  RIP college - HS growth gap (25-55): %.1f pp (paper: 43.4%%-42.7%%=0.7pp -> SAME)\n',...
    100*(results(3).cg55-results(4).cg55));
fprintf('  Data (CEX): college +74%%, HS +36%%\n');

%% Figure 11 analog: normalized mean consumption profiles, HIP (left) vs RIP (right)
fig=figure(1); clf;
subplot(1,2,1); hold on;
plot(25:65,results(1).meanc(1:Jwork+1)/results(1).meanc(1),'b-','LineWidth',1.5);
plot(25:65,results(2).meanc(1:Jwork+1)/results(2).meanc(1),'r--','LineWidth',1.5);
xlabel('age'); ylabel('consumption (age 25 = 1)'); title('HIP model'); legend('college','high school','Location','northwest');
subplot(1,2,2); hold on;
plot(25:65,results(3).meanc(1:Jwork+1)/results(3).meanc(1),'b-','LineWidth',1.5);
plot(25:65,results(4).meanc(1:Jwork+1)/results(4).meanc(1),'r--','LineWidth',1.5);
xlabel('age'); title('RIP model'); legend('college','high school','Location','northwest');
sgtitle('Figure 11 analog: avg life-cycle consumption by education (accurate solution)');
saveas(fig,'SavedOutput/Graphs/Replication_Fig11_education.png');
fprintf('Wrote SavedOutput/Graphs/Replication_Fig11_education.png\n');

save('RunEducationGroups_results.mat','results','edu','hipP','ripP','sab_rule','-v7.3');
fprintf('\n=== RunEducationGroups done ===\n');
diary off;
