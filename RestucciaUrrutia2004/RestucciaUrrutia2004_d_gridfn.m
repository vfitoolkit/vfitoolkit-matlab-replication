function d_grid=RestucciaUrrutia2004_d_gridfn(agej,g,gamma, sigma_b,rho_b, n_b,n_bhat, tauchenq,maxgridptparam)

% n_d=[n_bhatprime,1; n_s; 1]; % bhatprime,s

if agej==1
    [b_grid, ~]=TauchenMethod(0,sigma_b, rho_b, n_b, tauchenq);
    b_grid=exp(b_grid);
    minbhat=b_grid(1)*g^gamma;
    maxbhat=b_grid(end)*(maxgridptparam*g)^gamma;
    bhatprime_grid=exp(linspace(log(minbhat),log(maxbhat),n_bhat)');
    s_grid=[0;1];
elseif agej==2
    bhatprime_grid=1;
    s_grid=1;
end

d_grid=[bhatprime_grid; s_grid];

end