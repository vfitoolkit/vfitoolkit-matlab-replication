function StationaryDist=DiazGimenezPrescottAlvarezFitzgerald1992_StationaryDist(Policy,Params,n_d,n_a,n_s,n_z,n_sz,n_A_zero,pi_sz)

% Simulation of the model is highly non-standard due to the combination of
% stochastic death in infinitely lived agents together with their not
% caring about future generations. [Normally would either use
% finite-lifetime, or if using infinite-lived they would be dynasties
% that care about future generations.]

% Note that for n_z=1 this stationary distribution is also the stationary
% distribution in which all the agents actually find themselves. If n_z>1
% then it does not represent where all the agents actually at any point in
% time (which varies with z and must be simulated for a given path of z)
% but does still represent the asymptotic distribution of agents and so can
% be used to find asymptotic values of many model statistics of interest.

% The following makes use of some 'internal functions' of the VFI Toolkit
% to deal with the non-standard agent distribution. Most of it is simply a
% minor modification of contents of StationaryDist_Case1().
simoptions.parallel=2;
simoptions.tolerance=10^(-9);
simoptions.maxit=5*10^4;
PolicyKron=KronPolicyIndexes_Case1(Policy, n_d, n_a, n_sz,simoptions);

% Create a stationary dist. Am simply setting it up as if all newborns.
StationaryDist=zeros([n_a,n_sz],'gpuArray');
StationaryDist(n_A_zero,1,1,:)=Params.phi_1/prod(n_z);
StationaryDist(n_A_zero,1,2,:)=Params.phi_2/prod(n_z);

N_a=prod(n_a);
N_sz=prod(n_sz);
StationaryDistKron=reshape(StationaryDist,[N_a*N_sz,1]);

% simoptions.parallel==2 % Using the GPU
optaprime=reshape(PolicyKron(2,:,:),[1,N_a*N_sz]);

Ptemp=zeros(N_a,N_a*N_sz,'gpuArray');
Ptemp(optaprime+N_a*(gpuArray(0:1:N_a*N_sz-1)))=1;
Ptran=(kron(pi_sz',ones(N_a,N_a,'gpuArray'))).*(kron(ones(N_sz,1,'gpuArray'),Ptemp));

StationaryDistKronOld=zeros(N_a*N_sz,1,'gpuArray');
SScurrdist=sum(abs(StationaryDistKron-StationaryDistKronOld));
SScounter=0;

while SScurrdist>simoptions.tolerance && (100*SScounter)<simoptions.maxit
    
    for jj=1:100
        StationaryDistKron=Ptran*StationaryDistKron; %No point checking distance every single iteration. Do 100, then check.
        % NON-STANDARD PART: reallocate the dead to newborns
        StationaryDistKron=reshape(StationaryDistKron,[N_a,N_sz]); % Could use cleaver indexing, but feeling lazy, so just going to reshape and then un-reshape
        for z_c=1:n_z
            MassOfDead=sum(StationaryDistKron(:,n_s*z_c));
            StationaryDistKron(n_A_zero,1+n_s*(z_c-1))=StationaryDistKron(n_A_zero,1+n_s*(z_c-1))+Params.phi_1*MassOfDead; % want n_A_zero and K==1, but this is anyway just n_A_zero
            StationaryDistKron(n_A_zero,2+n_s*(z_c-1))=StationaryDistKron(n_A_zero,2+n_s*(z_c-1))+Params.phi_2*MassOfDead;
            StationaryDistKron(:,n_s*z_c)=0;
        end
        StationaryDistKron=reshape(StationaryDistKron,[N_a*N_sz,1]);
    end
    
    StationaryDistKronOld=StationaryDistKron;
    StationaryDistKron=Ptran*StationaryDistKron;
    % NON-STANDARD PART: reallocate the dead to newborns
    StationaryDistKron=reshape(StationaryDistKron,[N_a,N_sz]); % Could use cleaver indexing, but feeling lazy, so just going to reshape and then un-reshape
    for z_c=1:n_z
        MassOfDead=sum(StationaryDistKron(:,n_s*z_c));
        StationaryDistKron(n_A_zero,1+n_s*(z_c-1))=StationaryDistKron(n_A_zero,1+n_s*(z_c-1))+Params.phi_1*MassOfDead; % want n_A_zero and K==1, but this is anyway just n_A_zero
        StationaryDistKron(n_A_zero,2+n_s*(z_c-1))=StationaryDistKron(n_A_zero,2+n_s*(z_c-1))+Params.phi_2*MassOfDead;
        StationaryDistKron(:,n_s*z_c)=0;
    end
    StationaryDistKron=reshape(StationaryDistKron,[N_a*N_sz,1]);

    SScurrdist=sum(abs(StationaryDistKron-StationaryDistKronOld));    
    
    SScounter=SScounter+1;
    if rem(SScounter,50)==0
        SScounter
        SScurrdist
    end
end

StationaryDist=reshape(StationaryDistKron,[n_a,n_sz]);

sum(reshape(StationaryDist,[N_a,N_sz]),1)

end