% DFJ2010_Grids: asset grid (DFJ2010_AssetGrid), health and medical-expense shock grids, e-shock options, and the
% ptype-specific joint transition matrices pi_z_J. Script shared by DeNardiFrenchJones2010.m and DFJ2010_EndogMedex.m;
% needs n_z, n_e, N_j, Names_i, Params.rho_m, Params.sigma_mepsilon, Params.sigma_mxi, pi_h_J in the workspace.
% [Moved out of DeNardiFrenchJones2010.m on 2026-09-29, line for line.]

%% Grids
% a_grid and n_a are set by DFJ2010_AssetGrid (called with the state space above): 801 points on [0,300k] where all the
% target medians lie, 200 on (300k,2m]. [Until 2026-09-23: 2*(10^6)*linspace(0,1,501)'.^3, cubic on [0,2m].]
% note: borrowing constraint aprime>=0 is implicitly imposed by grid
% Note: DFJ2010 codes have max 'cash-on-hand' (assets +income) of
% 1.6million, 1million 200,000 depending on type. (lines 75-77 of their wealth20_Singles.c)

h_grid=[0;1;2;3]; % health status, 0=good, 1=bad, 2=death, 3=dead

% AR(1) markov medical expense shock, zeta
[zeta_grid,pi_zeta]=discretizeAR1_FarmerToda(0,Params.rho_m, Params.sigma_mepsilon, n_z(2));
% Note: DFJ2010 used Tauchen-Hussey (Guassian quadrature) to discretize
% [zeta_grid,pi_zeta]=discretizeAR1_TauchenHussey(0,Params.rho_m, Params.sigma_mepsilon, n_z(2)); % max is 1.26 instead of 1.73

% i.i.d. medical expense shock, xi
[xi_grid,pi_xi]=discretizeAR1_FarmerToda(0,0, Params.sigma_mxi, n_e);
% Note: DFJ2010 used Tauchen-Hussey (Guassian quadrature) to discretize
% [xi_grid,pi_xi]=discretizeAR1_TauchenHussey(0,0, Params.sigma_mxi, n_e); % max is 2.77 instead of 3.05
pi_xi=pi_xi(1,:)';

z_grid=[h_grid; zeta_grid];
e_grid=xi_grid; % [Until 2026-09-27 this line read e_grid=pi_xi, the probabilities: the iid shock xi then took values in (0,0.3) instead of the grid points on [-3.06,3.06], removing two thirds of the medical-expense shock variance. All results before that date are affected.]
pi_e=pi_xi;

d_grid=[];

% put e into vfoptions and simoptions
vfoptions.n_e=n_e;
vfoptions.e_grid=e_grid;
vfoptions.pi_e=pi_e;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;

% pi_z_J: combine pi_h_J and pi_zeta
% female
pi_z_J.femaleq1=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq1(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq1(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq2=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq2(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq2(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq3=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq3(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq3(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq4=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq4(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq4(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq5=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq5(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq5(:,:,jj)); % kron() in reverse order
end
% male
pi_z_J.maleq1=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq1(:,:,jj)=kron(pi_zeta,pi_h_J.maleq1(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq2=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq2(:,:,jj)=kron(pi_zeta,pi_h_J.maleq2(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq3=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq3(:,:,jj)=kron(pi_zeta,pi_h_J.maleq3(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq4=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq4(:,:,jj)=kron(pi_zeta,pi_h_J.maleq4(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq5=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq5(:,:,jj)=kron(pi_zeta,pi_h_J.maleq5(:,:,jj)); % kron() in reverse order
end

