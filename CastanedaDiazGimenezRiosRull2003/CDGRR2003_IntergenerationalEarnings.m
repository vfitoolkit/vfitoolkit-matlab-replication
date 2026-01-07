function IntergenerationalEarningsCorr=CDGRR2003_IntergenerationalEarnings(Nsim,StationaryDist, Policy, Params, simoptions, n_d,n_a,n_s,d_grid,a_grid,s_grid, pi_s, e,w,J)
% Intergenerational correlation of average lifetime earnings

% This statistic is calculated by simulation
% Simulate Nsim people
% At birth, people in the model are 21, we need to follow them till age 60,
% so follow them for T=40 periods

%Some variables we will later need
N_a=prod(n_a);
N_s=prod(n_s);

Policy=gather(Policy);
pi_s=gather(pi_s);

%%
N_z=N_s;
N_zprime=N_z;
% aprimeFnParamNames in same fashion
l_d2=length(n_d(2));
l_z=length(n_s);
temp=getAnonymousFnInputNames(simoptions.aprimeFn);
if length(temp)>(l_d2+2*l_z)
    aprimeFnParamNames={temp{l_d2+2*l_z+1:end}}; % the first inputs will always be (d2,a2)
else
    aprimeFnParamNames={};
end

whichisdforinheritasset=length(n_d);  % is just saying which is the decision variable that influences the experience asset (it is the 'last' decision variable)
aprimeFnParamsVec=CreateVectorFromParams(Params, aprimeFnParamNames);
[a2primeIndexes, a2primeProbs]=CreateaprimePolicyInheritanceAsset(Policy,simoptions.aprimeFn, whichisdforinheritasset, n_d, [],n_a, n_s, n_s, gpuArray(d_grid), a_grid, gpuArray(s_grid), gpuArray(s_grid), aprimeFnParamsVec);
% Note: aprimeIndexes and aprimeProbs are both [N_a,N_z,N_zprime]

Policy_a2prime=zeros(N_a,N_z,N_zprime,2,'gpuArray'); % the lower grid point
PolicyProbs=zeros(N_a,N_z,N_zprime,2,'gpuArray'); % The fourth dimension is lower/upper grid point
Policy_a2prime(:,:,:,1)=a2primeIndexes; % lower grid point
Policy_a2prime(:,:,:,2)=a2primeIndexes+1; % upper grid point
PolicyProbs(:,:,:,1)=a2primeProbs; % probability of lower grid point
PolicyProbs(:,:,:,2)=1-a2primeProbs; % probability of upper grid point
Policy_aprime=Policy_a2prime;

%%
% Policy_aprime and PolicyProbs are currently [N_a,N_z,N_zprime,N_probs]
Policy_aprimezprime=Policy_aprime+N_a*shiftdim(gpuArray(0:1:N_z-1),-1);  % Note: add z' index following the z' dimension
Policy_aprimezprime=gather(reshape(Policy_aprimezprime,[N_a*N_z,N_zprime*2])); % sparse() requires inputs to be 2-D
PolicyProbs=gather(reshape(PolicyProbs,[N_a*N_z,N_zprime*2])); % sparse() requires inputs to be 2-D

% Precompute
II2=repmat((1:1:N_a*N_z)',1,N_zprime*2); %  Index for this period (a,z), note the N_zprime*N_probs-copies

pi_s=sparse(gather(repelem(repmat(pi_s,1,2),N_a,1)));

% Transition matrix
P=sparse(Policy_aprimezprime,II2,PolicyProbs.*pi_s,N_a*N_z,N_a*N_z)'; % Note: sparse() will accumulate at repeated indices
% note the tranpose

% Now turn P into a cumulative distn
P=cumsum(P,2);

%%
StationaryDistKron=reshape(gather(StationaryDist),[N_a,N_s]);

% The actual calculations
AvgIndividualEarnings_Generation1=zeros(Nsim,1);
AvgIndividualEarnings_Generation2=zeros(Nsim,1);

CumSumSteadyStateDistVec_Old=cumsum(reshape(StationaryDistKron(:,J+1:2*J),[N_a*J,1]));
CumSumSteadyStateDistVec_Old=CumSumSteadyStateDistVec_Old./max(CumSumSteadyStateDistVec_Old); %renormalize so it is a cdf
parfor i=1:Nsim
    % First, t=1, they are born to a random position.
    % Observe that all old are equally likely to die, thus we use stationary
    % distribution to pick a 'random old position', and we then figure out
    % where this would leave the newborn (note, set from the old to newborn is random).
    [~,DyingIndex]=max(CumSumSteadyStateDistVec_Old>rand(1,1));
    temp=ind2sub_homemade([N_a,J],DyingIndex); a_c=temp(1); s_c_minus_J=temp(2);
    s_c=s_c_minus_J+J; %this is to compensate for the fact the DyingIndex is taken only from the old.
    % So our old person who is about to die is in [a_c,s_c].
    %currstate=sub2ind_homemade([N_a,N_s],[a_c,s_c]);
    % Thus our newborn will be born in exogenous state s_newborn (we have to 'force' death to occour):
    pi_temp=cumsum(pi_s(s_c,1:J));
    pi_temp=pi_temp./max(pi_temp);
    [~,s_newborn]=max(pi_temp>rand(1,1));
    % And with asset holdings a_newborn:
    a_newborn=a2primeIndexes(a_c,s_c,s_newborn);
    a_newborn=a_newborn+(a2primeProbs(a_c,s_c,s_newborn)<rand(1,1));
    
    % So newborn is in [a_newborn,s_newborn]
    AvgIndividualEarnings=e(s_newborn)*d_grid(Policy(1,a_newborn,s_newborn))*w;
    currstate=sub2ind_homemade([N_a,N_s],[a_newborn,s_newborn]);
    t=1;
    NumOfDeaths=0;
    while NumOfDeaths<2 %ie. for this generation and the next
        sold_c=s_c;
        t=t+1;
        [~,currstate]=max(P(currstate,:)>rand(1,1));
        temp=ind2sub_homemade([N_a,N_s],currstate); a_c=temp(1); s_c=temp(2);
        if J+1<=sold_c<=2*J && s_c<=J %if individual dies
            NumOfDeaths=NumOfDeaths+1;
            if NumOfDeaths==1
                AvgIndividualEarnings_Generation1(i)=AvgIndividualEarnings;
            else
                AvgIndividualEarnings_Generation2(i)=AvgIndividualEarnings;
            end
            AvgIndividualEarnings=e(s_c)*d_grid(Policy(1,a_c,s_c))*w;
            t=1;
        elseif s_c<=J %If they are still working update their earnings
            AvgIndividualEarnings=((t-1)/t)*AvgIndividualEarnings+(1/t)*e(s_c)*d_grid(Policy(1,a_c,s_c))*w;
        end
    end
end
AvgIndividualEarnings_Generations=[AvgIndividualEarnings_Generation1, AvgIndividualEarnings_Generation2];

Cov_old_young=sum(AvgIndividualEarnings_Generations(:,1).*AvgIndividualEarnings_Generations(:,2))/Nsim-(sum(AvgIndividualEarnings_Generations(:,1))/Nsim)*(sum(AvgIndividualEarnings_Generations(:,2))/Nsim);
StdDev_old=sqrt(sum(AvgIndividualEarnings_Generations(:,2).^2)/Nsim-(sum(AvgIndividualEarnings_Generations(:,2))/Nsim)^2);
StdDev_young=sqrt(sum(AvgIndividualEarnings_Generations(:,1).^2)/Nsim-(sum(AvgIndividualEarnings_Generations(:,1))/Nsim)^2);
IntergenerationalEarningsCorr=Cov_old_young/(StdDev_old*StdDev_young);

end