% DFJ2010_InitialDistCohorts: the 1996 initial distribution by ptype and cohort from SavedOutput/DFJ2010_DataWork.mat,
% cohort entry ages, jequaloneDist with an N_j dimension, age weights, ptype masses, and the aliveq1..5 restrictions.
% Script shared by DeNardiFrenchJones2010.m and DFJ2010_EndogMedex.m; needs n_a, a_grid, n_z, n_e, N_j, Names_i,
% Params.agejshifter, pi_zeta, simoptions in the workspace. [Moved out of DeNardiFrenchJones2010.m on 2026-09-29.]

%% Setup initial agent distribution: the 1996 AHEAD data, by ptype and by cohort (from DFJ2010_DataWork.m)
% DFJ2010 start each simulated individual from the 1996 data distribution of (age, sex, PI, health, assets).
% DFJ2010_DataWork.m puts that into InitialDist.InitDist.(ptype), an n_a x 2 x 33 mass over (asset grid point,
% health {good,bad}, 1996 age 70..102), and InitialDist.ptypemasses (1996 shares of the ten ptypes).
% The toolkit handles cohorts entering at different ages through jequaloneDist with an N_j dimension, together with
% estimoptions.cohortagejshifter (estimation) or simoptions.jequaloneDistAge (one cohort at a time).
% DFJ2010 group the 1996 sample into five birth cohorts (1996 ages 72-76, 77-81, 82-86, 87-91, 92-102) and match each
% cohort's median assets in 1998,...,2006. Here each cohort enters the model at its rounded average 1996 age (DFJ2010's
% own avgage96), with its members' assets and health pooled over the cohort's ages.
if ~exist('./SavedOutput/DFJ2010_DataWork.mat','file')
    error('Run DFJ2010_DataWork.m first: it builds SavedOutput/DFJ2010_DataWork.mat with the 1996 initial distribution and the GMM targets')
end
DataWork=load('./SavedOutput/DFJ2010_DataWork.mat','Targets','InitialDist');
if DataWork.InitialDist.n_a~=n_a || max(abs(DataWork.InitialDist.a_grid-a_grid))>0
    error('The asset grid in SavedOutput/DFJ2010_DataWork.mat differs from the model grid (DFJ2010_AssetGrid): rerun the data work, doPart(1)=1')
end
% Cohorts are defined on 1996 age (GAUSS getchrt on 1994 age = 1996 age - 2, with edges cohortedges=[69,74,79,84,89,102])
cohortedges=DataWork.Targets.cohortedges;
ncohorts=5;
age1996=70:102; % the 33 age slices of InitDist
cohortofage=zeros(1,N_j); % cohort of each 1996 age (0 for ages 70-71, which are in no DFJ cohort)
for c=1:ncohorts
    cohortofage((age1996-2)>cohortedges(c) & (age1996-2)<=cohortedges(c+1))=c;
end
% Entry age of each cohort: rounded average 1996 age of its members, all ptypes together
InitObs=DataWork.InitialDist.InitObs;
cohortentryage=zeros(1,ncohorts);
for c=1:ncohorts
    cohortentryage(c)=round(mean(InitObs.age1996(cohortofage(InitObs.age1996-Params.agejshifter)==c)));
end
fprintf('Cohorts enter the model at ages %s (DFJ2010: 74, 79, 84, 89, and about 95) \n',mat2str(cohortentryage))
% Medical expense shocks in 1996: DFJ2010 draw zeta from its stationary distribution and xi iid (GAUSS initmx)
zetastat=pi_zeta^1000; zetastat=zetastat(1,:)';
% jequaloneDist by ptype: [n_a, n_z(1), n_z(2), n_e, N_j], mass one per ptype. The slice at a cohort's entry age is that
% cohort's (assets, health) distribution from the 1996 data (health 0/1 are the first two of the four h states).
jequaloneDist=struct();
for ii=1:length(Names_i)
    D=DataWork.InitialDist.InitDist.(Names_i{ii}); % n_a x 2 x 33
    temp=zeros([n_a,n_z(1),n_z(2),n_e,N_j]);
    for c=1:ncohorts
        Dc=sum(D(:,:,cohortofage==c),3); % n_a x 2: mass of this cohort within this ptype
        temp(:,1:2,:,:,cohortentryage(c)-Params.agejshifter)=reshape(Dc,[n_a,2,1,1]).*reshape(zetastat,[1,1,n_z(2),1]).*reshape(pi_xi,[1,1,1,n_e]);
    end
    jequaloneDist.(Names_i{ii})=gpuArray(temp/sum(temp(:)));
end
clear temp D Dc

% Age weights: within-age statistics do not depend on them; the cohort machinery gives each cohort its entry-age weight at every age
AgeWeightsParamNames={'mewj'};
Params.mewj=ones(1,N_j)/N_j;

% Permanent type masses: the 1996 shares of the ten ptypes (women are 85 percent of the sample)
PTypeDistParamNames={'ptypemasses'};
Params.ptypemasses=DataWork.InitialDist.ptypemasses;
% Permanent-income quintile of each ptype, used by the aliveq1..aliveq5 conditional restrictions
% (the DFJ2010 targets pool men and women within a quintile)
Params.quintile.femaleq1=1; Params.quintile.femaleq2=2; Params.quintile.femaleq3=3; Params.quintile.femaleq4=4; Params.quintile.femaleq5=5;
Params.quintile.maleq1=1;   Params.quintile.maleq2=2;   Params.quintile.maleq3=3;   Params.quintile.maleq4=4;   Params.quintile.maleq5=5;
% alive and in permanent-income quintile q (pools the female and male ptypes of that quintile; these are the GMM targets)
simoptions.conditionalrestrictions.aliveq1=@(aprime,a,h,zeta,xi,quintile) ((h==0) || (h==1)) && (quintile==1);
simoptions.conditionalrestrictions.aliveq2=@(aprime,a,h,zeta,xi,quintile) ((h==0) || (h==1)) && (quintile==2);
simoptions.conditionalrestrictions.aliveq3=@(aprime,a,h,zeta,xi,quintile) ((h==0) || (h==1)) && (quintile==3);
simoptions.conditionalrestrictions.aliveq4=@(aprime,a,h,zeta,xi,quintile) ((h==0) || (h==1)) && (quintile==4);
simoptions.conditionalrestrictions.aliveq5=@(aprime,a,h,zeta,xi,quintile) ((h==0) || (h==1)) && (quintile==5);


