function IntergenerationalEarningsCorr=CDGRR2003_IntergenerationalEarnings(Nsim,StationaryDist, Policy, Phi_aprimeKron, n_d,n_a,n_s,d_grid, pi_s, e,w,J)
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
Phi_aprimeKron=gather(Phi_aprimeKron);

P=zeros(N_a,N_s,N_a,N_s); %P(a,z,aprime,zprime)=proby of going to (a',z') given in (a,z)
for a_c=1:N_a
    for s_c=1:N_s
        for sprime_c=1:N_s
            d2_c=Policy(1,a_c,s_c);
            optaprime=Phi_aprimeKron(d2_c,sprime_c,s_c);
            P(a_c,s_c,optaprime,sprime_c)=gather(pi_s(s_c,sprime_c)/sum(pi_s(s_c,:)));
        end
    end
end
P=reshape(P,[N_a*N_s,N_a*N_s]);
%Now turn P into a cumulative distn
P=cumsum(P,2);

StationaryDistKron=reshape(gather(StationaryDist),[N_a,N_s]);

% The actual calculations
AvgIndividualEarnings_Generation1=zeros(Nsim,1);
AvgIndividualEarnings_Generation2=zeros(Nsim,1);

CumSumSteadyStateDistVec_Old=cumsum(reshape(StationaryDistKron(:,J+1:2*J),[N_a*J,1]));
CumSumSteadyStateDistVec_Old=CumSumSteadyStateDistVec_Old./max(CumSumSteadyStateDistVec_Old); %renormalize so it is a cdf
parfor i=1:Nsim
    %First, t=1, they are born to a random position.
    %Observe that all old are equally likely to die, thus we use stationary
    %distribution to pick a 'random old position', and we then figure out
    %where this would leave the newborn (note, set from the old to newborn is random).
    [~,DyingIndex]=max(CumSumSteadyStateDistVec_Old>rand(1,1));
    temp=ind2sub_homemade([N_a,J],DyingIndex); a_c=temp(1); s_c_minus_J=temp(2);
    s_c=s_c_minus_J+J; %this is to compensate for the fact the DyingIndex is taken only from the old.
    %So our old person who is about to die is in [a_c,s_c].
    %currstate=sub2ind_homemade([N_a,N_s],[a_c,s_c]);
    %Thus our newborn will be born in exogenous state s_newborn (we have to 'force' death to occour):
    pi_temp=cumsum(pi_s(s_c,1:J));
    pi_temp=pi_temp./max(pi_temp);
    [~,s_newborn]=max(pi_temp>rand(1,1));
    %And with asset holdings a_newborn:
    d2_c=Policy(2,a_c,s_c);
    a_newborn=Phi_aprimeKron(d2_c,s_newborn,s_c);
    %So newborn is in [a_newborn,s_newborn]
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