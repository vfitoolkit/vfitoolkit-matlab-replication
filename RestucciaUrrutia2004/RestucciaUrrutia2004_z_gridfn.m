function [z_grid,pi_z]=RestucciaUrrutia2004_z_gridfn(agej, sigma_b,rho_b,earningsshocksize, n_b,n_thetahat,n_earnings, tauchenq)

% The actual z_grid on 'innate ability' is same for all ages, but transitions are only possible between generations. 
if agej==1
    [b_grid, ~]=TauchenMethod(0,(1-rho_b^2)*sigma_b^2, rho_b, n_b, tauchenq);
    b_grid=exp(b_grid);
    pi_b=eye(n_b,n_b,'gpuArray');
    thetahat_grid=1;
    pi_thetahat=ones(1,n_thetahat)/n_thetahat; % Uniform distribution
    if n_earnings==1
        earnings_grid=1;
        pi_earnings=1;
    elseif n_earnings==2
        earnings_grid=[1-earningsshocksize; 1+earningsshocksize];
        pi_earnings=[1,0;0,1];
    end
elseif agej==2
    [b_grid, pi_b]=TauchenMethod(0,(1-rho_b^2)*sigma_b^2, rho_b, n_b, tauchenq);
    b_grid=exp(b_grid);
    thetahat_grid=linspace(0,1-1/n_thetahat,n_thetahat)'; % I deliberately make the top point less than one, so q(bhat) above this have 'probability one of college completion'. (round the 99%ers up rather than down)
    pi_thetahat=ones(n_thetahat,1); % Goes to the only next period value with certainty
    if n_earnings==1
        earnings_grid=1;
        pi_earnings=1;
    elseif n_earnings==2
        earnings_grid=[1-earningsshocksize; 1+earningsshocksize];
        pi_earnings=[0.5,0.5;0.5,0.5];
    end
end
% I have not enforced the mean(b)=1, so this might be slightly out.

z_grid=[b_grid; thetahat_grid; earnings_grid];
pi_z=kron(pi_earnings, kron(pi_thetahat,pi_b)); % Note, this always is done in 'reverse order'

end