%% Takes a closer look at the inc_process.mat from the Guerrieri & Lorenzoni (2017) codes.
% This is how I figured out that they use tauchenq=2.1, and that this is
% based on targeting that sigmasqz takes the correct value.

%% Load the log-income process (which is 12-state Tauchen method approximation of AR(1))
load ./PaperMaterials/replication-codes-for-Credit-Crises-2017-1f7cb32/inc_process.mat

% Simulate a time series
burnin=10^3;
T=10^4
currenttimeseries_index=1
cumPr=cumsum(Pr,2);
% Start with the indexes
for ii=1:burnin
    [~,currenttimeseries_index]=max(cumPr(currenttimeseries_index,:)>rand(1,1));
end
for ii=1:T
    [~,currenttimeseries_index]=max(cumPr(currenttimeseries_index,:)>rand(1,1));
    timeseries_index(ii)=currenttimeseries_index;
end
% Now the values
timeseries_values=x(timeseries_index);
% Variance of these
var(timeseries_values)

% It should be
rho=0.967;
sigmasq_epsilon=0.017;
sigmasq_z=sigmasq_epsilon/(1-rho^2)


%% Use Tauchen Method
tauchenoptions.parallel=1;
Params.tauchenq=2.1
n_theta=13;
Params.rho=rho;
Params.sigmasq_epsilon=sigmasq_epsilon;
[theta1_grid, pi_theta1]=TauchenMethod(0, Params.sigmasq_epsilon, Params.rho, n_theta-1, Params.tauchenq,tauchenoptions);

cumpi_theta1=cumsum(pi_theta1,2);
% Start with the indexes
for ii=1:burnin
    [~,currenttimeseries_index]=max(cumpi_theta1(currenttimeseries_index,:)>rand(1,1));
end
for ii=1:T
    [~,currenttimeseries_index]=max(cumpi_theta1(currenttimeseries_index,:)>rand(1,1));
    timeseries_index(ii)=currenttimeseries_index;
end
% Now the values
timeseries_values=theta1_grid(timeseries_index);
% Variance of these
var(timeseries_values)
