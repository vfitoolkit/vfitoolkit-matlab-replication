% DFJ2010_EndogMedex4: the endogenous medical-expense model solved with medical goods as a DECISION VARIABLE on a grid
% (n_d=101, d=m), to compare runtime, GPU memory and the solution with the bisection approach of DFJ2010EM_ReturnFn
% (part 6). Called from DeNardiFrenchJones2010.m as doPart(8); needs part 6 in the same run (the endogenous parameters,
% grids, the bisection V and Policy, the cohorts and the targets).
%
% The grid for m is log-spaced from $10 to $2 million (101 points, about 13 percent apart): with omega>1, m=0 gives
% utility -Inf so it is excluded. The comparison is: value function (the grid solution can only be weakly worse),
% the asset policy, and the five-cohort asset-median objective and mean out-of-pocket medical expenses by cohort x
% PI half x year, which are what the estimation uses.
fprintf('\n===== Part 8: medical goods as a decision variable (n_d=101) vs the bisection split =====\n')
n_d_EM=101;
d_grid_EM=exp(linspace(log(10),log(2e6),n_d_EM))';
ReturnFn_dgrid=@(d,aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    DFJ2010EM_ReturnFn_dgrid(d,aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6);
vfoptions_d=vfoptions;
vfoptions_d.verbose=0;
vfoptions_d.verboseparams=0;
gpu=gpuDevice;
fprintf('GPU memory available before the solve: %.2f GB (of %.2f GB)\n',gpu.AvailableMemory/1e9,gpu.TotalMemory/1e9)
% Try the default memory setting first; if the return matrix (n_d x n_a x level1n x 420 shock states) does not fit, loop over e
vfoptions_d.lowmemory=0;
tic;
try
    [V_d,Policy_d]=ValueFnIter_Case1_FHorz_PType(n_d_EM,n_a,n_z,N_j,Names_i,d_grid_EM,a_grid,z_grid,pi_z_J,ReturnFn_dgrid,Params,DiscountFactorParamNames,vfoptions_d);
    lowmemused=0;
catch ME
    fprintf('Solve with vfoptions.lowmemory=0 failed (%s); retrying with lowmemory=1 (loop over e)\n',ME.message)
    vfoptions_d.lowmemory=1;
    tic;
    [V_d,Policy_d]=ValueFnIter_Case1_FHorz_PType(n_d_EM,n_a,n_z,N_j,Names_i,d_grid_EM,a_grid,z_grid,pi_z_J,ReturnFn_dgrid,Params,DiscountFactorParamNames,vfoptions_d);
    lowmemused=1;
end
vftime_d=toc;
gpu=gpuDevice;
fprintf('Value fn with m on a %d-point grid (lowmemory=%d) took %.1f seconds; bisection split took %.1f seconds; GPU memory available after: %.2f GB\n',n_d_EM,lowmemused,vftime_d,vftime,gpu.AvailableMemory/1e9)

% Value function and asset policy comparison, alive states only (ptype by ptype)
fprintf('Value fn and asset policy, grid-m vs bisection, over the alive states (h=0,1) at all ages:\n')
% asset policy values for all types at once: the bisection Policy has the aprime index (+ L2 rows), the grid-m Policy has the d index first
PVb=PolicyInd2Val_FHorz_PType(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,simoptions);
PVg=PolicyInd2Val_FHorz_PType(Policy_d,n_d_EM,n_a,n_z,N_j,d_grid_EM,a_grid,simoptions);
for ii=1:length(Names_i)
    Vb=V.(Names_i{ii}); Vg=V_d.(Names_i{ii}); % [n_a,4,15,7,N_j]
    alive=false(size(Vb)); alive(:,1:2,:,:,:)=true;
    dV=gather(Vg(alive)-Vb(alive)); Vbb=gather(Vb(alive));
    apb=gather(PVb.(Names_i{ii})(1,:,:,:,:,:)); apg=gather(PVg.(Names_i{ii})(2,:,:,:,:,:)); % aprime is row 1 (no d) / row 2 (after d)
    apb=apb(:); apg=apg(:);
    fprintf('   %-9s V: max |diff| %.3e, mean |diff|/|V| %.3e, share of states where grid V is lower %.3f; aprime: max |diff| %.0f dollars, mean |diff| %.0f\n',Names_i{ii},max(abs(dV)),mean(abs(dV)./abs(Vbb)),mean(dV<-1e-12),max(abs(apg-apb)),mean(abs(apg-apb)))
end

% The five-cohort statistics the estimation uses, with the grid-m policy
FnsToEvaluate_d.assets=@(d,aprime,a,h,zeta,xi) a;
FnsToEvaluate_d.oopmedex=@(d,aprime,a,h,zeta,xi,q_healthbad,q_healthgood) max(((h==0)*q_healthgood+(h==1)*q_healthbad)*d,250); % gross out-of-pocket on the grid, bottom-coded as the data
simoptions_d=simoptions_val;
modelmedian_d=nan(Jdata,1);
ModelMedex_d.Mean=nan(5,2,6);
Params_c=Params;
tic;
for c=1:ncohorts
    jequaloneDist_c=struct();
    cohortmass_c=zeros(1,length(Names_i));
    for ii=1:length(Names_i)
        temp=jequaloneDist.(Names_i{ii})(:,:,:,:,cohortentryage(c)-Params.agejshifter);
        cohortmass_c(ii)=gather(sum(temp(:)));
        if cohortmass_c(ii)>0
            jequaloneDist_c.(Names_i{ii})=temp/cohortmass_c(ii);
        else
            jequaloneDist_c.(Names_i{ii})=ones(size(temp),'gpuArray')/numel(temp);
        end
    end
    Params_c.ptypemasses=Params.ptypemasses(:)'.*cohortmass_c;
    Params_c.ptypemasses=Params_c.ptypemasses/sum(Params_c.ptypemasses);
    Params_c.mewj=Params.mewj(cohortentryage(c)-Params.agejshifter)*ones(1,N_j);
    simoptions_d.jequaloneDistAge=cohortentryage(c)-Params.agejshifter;
    StationaryDist_c=StationaryDist_Case1_FHorz_PType(jequaloneDist_c,AgeWeightsParamNames,PTypeDistParamNames,Policy_d,n_d_EM,n_a,n_z,N_j,Names_i,pi_z_J,Params_c,simoptions_d);
    ACS_c=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_c,Policy_d,FnsToEvaluate_d,Params_c,n_d_EM,n_a,n_z,N_j,Names_i,d_grid_EM,a_grid,z_grid,simoptions_d);
    rows=find(MomentList.cohort==c);
    for rr=rows'
        modelmedian_d(rr)=gather(ACS_c.(['aliveq',num2str(MomentList.PIq(rr))]).assets.Median(cohortentryage(c)+MomentList.year(rr)-1996-Params.agejshifter));
    end
    for hh=1:2
        for y=2:6
            jj=cohortentryage(c)+years(y)-1996-Params.agejshifter;
            if jj<=N_j
                ModelMedex_d.Mean(c,hh,y)=gather(ACS_c.(['aliveh',num2str(hh)]).oopmedex.Mean(jj));
            end
        end
    end
end
cohortlooptime_d=toc;
gap_d=modelmedian_d(tkorder)-MomentList.median(tkorder);
fprintf('Five-cohort loop with the grid-m policy took %.1f seconds (bisection: %.1f)\n',cohortlooptime_d,cohortlooptime)
fprintf('Asset-median objective: grid-m %.1f vs bisection %.1f; max |difference| in the 102 model medians %.0f dollars\n',gap_d'*WeightingMatrix*gap_d,gap'*WeightingMatrix*gap,max(abs(modelmedian_d-modelmedian)))
okm=~isnan(ModelMedex.Mean(:,:,2:6)) & ~isnan(ModelMedex_d.Mean(:,:,2:6));
mb=ModelMedex.Mean(:,:,2:6); mg=ModelMedex_d.Mean(:,:,2:6);
fprintf('Mean out-of-pocket medex over the cells: bisection %.0f vs grid-m %.0f; max |difference| %.0f dollars (the grid-m one is gross of Medicaid and on a 13-percent grid)\n',mean(mb(okm)),mean(mg(okm)),max(abs(mb(okm)-mg(okm))))
save ./SavedOutput/DFJ2010_EndogMedex_dgrid.mat vftime_d vftime lowmemused modelmedian_d modelmedian ModelMedex_d ModelMedex cohortlooptime_d cohortlooptime n_d_EM d_grid_EM
