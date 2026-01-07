function RatioEarningsOldYoung=CDGRR2003_RatioEarningsOldYoung(Nsim, StationaryDist, Policy, Params, simoptions, n_d,n_a,n_s, d_grid,a_grid,s_grid, pi_s, e,w,J)

% Life-cycle profile of earnings
%  (avg. earnings of HH aged 46 to 50)/(avg. earnings of HH aged 26 to 30)

% This statistic is calculated by simulation
% Simulate Nsim people
% At birth, people in the model are 21, we follow them till age 50,
% so follow them for T=30 periods
% Only people who work for all 40 periods are considered in this sample.
% Anyone who retires at any point is deleted from sample and replaced.

K=3; % I am going to simulate 'K' times as many people as needed, and then just keep all of those who actually did retire

% Some variables we will later need
N_a=prod(n_a);
N_s=prod(n_s);


StationaryDistKron=reshape(StationaryDist,[N_a,N_s]);
CumSumSteadyStateDistVec_Old=cumsum(reshape(StationaryDistKron(:,J+1:2*J),[N_a*J,1]));
CumSumSteadyStateDistVec_Old=gather(CumSumSteadyStateDistVec_Old./max(CumSumSteadyStateDistVec_Old)); %renormalize so it is a cdf
Policy=gather(Policy);
pi_s=gather(pi_s);


% The actual calculations
% Note that the sample is taken such that only people who remain workers
% for 40 years are considered (if they retire they are deleted from this
% sample).
IndividualEarnings_Young=zeros(K*Nsim,1); %age 26-30
IndividualEarnings_Old=zeros(K*Nsim,1);   %age 46-50
cumsumpi_s=cumsum(pi_s(1:J,:),2); %rows have to add to one anyway

pi_temp=cumsum(pi_s(J+1:2*J,1:J),2);
pi_temp=pi_temp./(max(pi_temp,[],2)*ones(1,J));

%%
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


%%
parfor i=1:K*Nsim 
    TheyRetired=0;
    YoungEarnings=zeros(1,20);
    OldEarnings=zeros(1,10);
    % First, t=1, they are born to a random position.
    % Observe that all old are equally likely to die, thus we use stationary distribution to pick a 'random old position', and we then figure out where this would leave the newborn (note, set from the old to newborn is random).
    [~,DyingIndex]=max(CumSumSteadyStateDistVec_Old>rand(1,1));
    temp=ind2sub_homemade([N_a,J],DyingIndex); a_c=temp(1); s_c_minus_J=temp(2);
    s_c=s_c_minus_J+J; % this is to compensate for the fact the DyingIndex is taken only from the old.
    % So our old person who is about to die is in [a_c,s_c].

    % Thus our newborn will be born in exogenous state s_newborn (we have to 'force' death to occour):
    [~,s_newborn]=max(pi_temp(s_c_minus_J,:)>rand(1,1));
    % And with asset holdings a_newborn:
    a_newborn=a2primeIndexes(a_c,s_c,s_newborn);
    a_newborn=a_newborn+(a2primeProbs(a_c,s_c,s_newborn)<rand(1,1));
        
    % So newborn is in [a_newborn,s_newborn]
    YoungEarnings(1)=e(s_newborn)*d_grid(Policy(1,a_newborn,s_newborn))*w;
    a_c=a_newborn;
    s_c=s_newborn;
    for t=2:20
        if s_c>J %Check that they are still working
            TheyRetired=1; %If not working mark them as having retired.
        else
            d2_c=Policy(2,a_c,s_c);
            a_c=d2_c; %a_c=Phi_aprimeKron(d1_c,s_c,new_s_c); %Actually, can ignore the estate tax, since we are only looking at people who don't retire (and so cannot die) anyway
            [~,s_c]=max(cumsumpi_s(s_c,:)>rand(1,1));
            YoungEarnings(t)=e(s_c)*d_grid(Policy(1,a_c,s_c))*w;
        end
        
    end
    if TheyRetired==0
        for t=21:30
            
            if s_c>J %Check that they are still working
                TheyRetired=1; %If not working mark them as having retired.
            else
                d2_c=Policy(2,a_c,s_c);
                a_c=d2_c; %a_c=Phi_aprimeKron(d1_c,s_c,new_s_c); %Actually, can ignore the estate tax, since we are only looking at people who don't retire (and so cannot die) anyway
                [~,s_c]=max(cumsumpi_s(s_c,:)>rand(1,1));
                OldEarnings(t-20)=e(s_c)*d_grid(Policy(1,a_c,s_c))*w;
            end
        end
    end
    if TheyRetired==0
        IndividualEarnings_Young(i)=mean(YoungEarnings(6:10)); % YoungEarnings(6:10) corresponds to ages 26-30
        IndividualEarnings_Old(i)=mean(OldEarnings(6:10)); % OldEarnings(6:10) corresponds to ages 46-50
%         disp('Got one!')
    end
end
IndividualEarnings_Young=IndividualEarnings_Young(IndividualEarnings_Young~=0);
IndividualEarnings_Old=IndividualEarnings_Old(IndividualEarnings_Old~=0);
disp('Number of sample points for RatioEarningsOldYoung')
size(IndividualEarnings_Old)
RatioEarningsOldYoung=mean(IndividualEarnings_Old)/mean(IndividualEarnings_Young);

%One oddity here is that since aging is random, it is possible that 21 to
%60 year olds may be retired. We do not wish to allow for this when
%calculating this statistic. So we only consider people who make it through
%the 40 periods without retiring.

end