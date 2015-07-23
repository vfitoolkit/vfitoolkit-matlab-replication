% Replicate Tauchen (1986) - Finite State Markov-Chain Approximations to Univariate and Vector Auotregressions
% REPLICATION IS INCOMPLETE. NEED TO GENERATE THE TABLES AND DO THE REGRESSIONS

tauchenoptions.parallel=0;

% Demonstrates how to use TauchenMethod and TauchenMethodVAR functions of VFI Toolkit


%% The AR(1) approximation
znum=9 % Tauchen uses 9 % Number of points must be odd (also tries out 5; see his Table 1)
T=10^5 % Tauchen uses 50
q=3 % Tauchen uses 3

% Use command TauchenMethod(mew,sigmasq,rho,znum,q,tauchenoptions)
mew=0;
rho_vec=[0.1,0.8,0.9]; % See Tauchen Table 1
sigmasq_vec=[0.101,0.167,0.229].^2; % See Tauchen Table 1

Table1Case=1; % 1,2, or 3
rho=rho_vec(Table1Case);
sigmasq=sigmasq_vec(Table1Case);

[z_grid, pi_z]=TauchenMethod(mew,sigmasq,rho,znum,q,tauchenoptions);


% Simulate the resulting markov chain
burnin=1000;
cumsum_pi_z=cumsum(pi_z,2);
z_simindexes=zeros(T,1);
z_simindexes(1)=ceil(znum/2);
for t=2:burnin
    temp_cumsum_pi_z=cumsum_pi_z(z_simindexes(1),:);
    temp_cumsum_pi_z(temp_cumsum_pi_z<=rand(1,1))=2;
    [~,z_simindexes(1)]=min(temp_cumsum_pi_z);
end
for t=2:T
    temp_cumsum_pi_z=cumsum_pi_z(z_simindexes(t-1),:);
    temp_cumsum_pi_z(temp_cumsum_pi_z<=rand(1,1))=2;
    [~,z_simindexes(t)]=min(temp_cumsum_pi_z);
end
z_simvals=nan(1,T);
for t=1:T
    z_simvals(t)=z_grid(z_simindexes(t));
end

% Check the mean values of the simulation
mean(z_simvals)
% They should be equal to the analytical
((1-rho)^(-1))*mew

%% The VAR(1) approximation
znum=[25;25] % Tauchen uses [9;9] % Number of points in each dimension must be odd.
T=10^5 % Tauchen uses 50
q=[3;3] % Tauchen uses [3;3]

% Use command TauchenMethodVAR(mew,sigmasq,rho,znum,q,tauchenoptions)
% Create states vector and transition matrix for the discrete markov process approximation of n-variable VAR(1) process z'=mew+rho*z+e, e~N(0,sigmasq), by Tauchens method
% (Abuse of notation: sigmasq in codes is a column vector, in the VAR notation it is a matrix in which all of the non-diagonal elements are zeros) 
mew=[0;0];
rho=[0.7,0.3; 0.2,0.5];
sigmasq=[0.1;0.1];

[z_grid, pi_z]=TauchenMethodVAR(mew,sigmasq,rho,znum,q,tauchenoptions);

% Simulate the resulting markov chain
burnin=1000;
cumsum_pi_z=cumsum(pi_z,2);
z_simindexes=zeros(T,1);
z_simindexes(1)=ceil(prod(znum)/2);
for t=2:burnin
    temp_cumsum_pi_z=cumsum_pi_z(z_simindexes(1),:);
    temp_cumsum_pi_z(temp_cumsum_pi_z<=rand(1,1))=2;
    [~,z_simindexes(1)]=min(temp_cumsum_pi_z);
end
for t=2:T
    temp_cumsum_pi_z=cumsum_pi_z(z_simindexes(t-1),:);
    temp_cumsum_pi_z(temp_cumsum_pi_z<=rand(1,1))=2;
    [~,z_simindexes(t)]=min(temp_cumsum_pi_z);
end
z_simvals=nan(2,T);
for t=1:T
    ind=z_simindexes(t);
    sub=ind2sub_homemade([znum(1),znum(2)],ind);
    z_simvals(1,t)=z_grid(sub(1));
    z_simvals(2,t)=z_grid(sub(2));
end

% Check the mean values of the simulation
mean(z_simvals')
% They should be equal to the analytical
((eye(2,2)-rho)^(-1))*mew

