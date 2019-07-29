function a_grid=RestucciaUrrutia2004_a_gridfn(agej,g,gamma,pupperbar, sigma_b,rho_b, n_b,n_bhat,n_h, tauchenq, maxgridptparam)

% n_a=[n_h, n_h; 1, n_bhat;1,n_s]; % h,bhat,s

% I use 8g rather than 5g used by RU2004 to set top of bhat as I had people
% running into the top of the grid on h (have called this maxgridptparam)

if agej==1
    [b_grid, ~]=TauchenMethod(0,sigma_b, rho_b, n_b, tauchenq);
    b_grid=exp(b_grid);
    minh=b_grid(1)*g^gamma;
    maxh=pupperbar*b_grid(end)*(maxgridptparam*g)^gamma;
    h_grid=exp(linspace(log(minh),log(maxh),n_h)');
    bhat_grid=1;
    s_grid=1;
elseif agej==2
    [b_grid, ~]=TauchenMethod(0,sigma_b, rho_b, n_b, tauchenq);
    b_grid=exp(b_grid);
    minbhat=b_grid(1)*g^gamma;
    maxbhat=b_grid(end)*(maxgridptparam*g)^gamma;
    minh=b_grid(1)*g^gamma; % =minbhat
    maxh=pupperbar*b_grid(end)*(maxgridptparam*g)^gamma; % = pupperbar*maxbhat;
    h_grid=exp(linspace(log(minh),log(maxh),n_h)');
    bhat_grid=exp(linspace(log(minbhat),log(maxbhat),n_bhat)');
    s_grid=[0;1];
end

a_grid=[h_grid; bhat_grid; s_grid];

end