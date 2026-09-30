function [k,n,output] = B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar)
% B2021_EntrepreneurStaticProblem finds optimal entrepreneurial choices given a collateral constraint.
%
% INPUTS:
%   a       : Assets (state variable)
%   theta   : Entrepreneurial ability (state variable)
%   r, w    : Interest rate and wage
%   lambda  : Collateral constraint parameter (k <= lambda * a)
%   delta   : Depreciation rate
%   gamma   : Capital share in production
%   upsilon : Span-of-control parameter
%   lbar    : Fixed own-labor input supplied by entrepreneurs
%
% OUTPUTS:
%   k      : Optimal capital demand
%   n      : Optimal hired labor demand
%   output : Optimal output
%
% MODEL:
%   Production function:
%       y = theta * (k^gamma * nbar^(1 - gamma))^upsilon
%
%   Profit maximized:
%       y - (r + delta) * k - w * n
%   subject to:
%       k <= lambda * a
%       n >= 0
%
% The function first solves the candidate interior solution with n > 0 and
% then re-solves under n = 0 whenever hired labor would otherwise be
% negative.

% Case 1: Unconstrained solution with n > 0

% Unconstrained level of capital (n > 0 candidate)
aux1   = (gamma / (r + delta))^(1 + (gamma * upsilon) / (1 - upsilon));
aux2   = ((1 - gamma) / w)^((upsilon * (1 - gamma)) / (1 - upsilon));
k_unc  = (theta * upsilon)^(1 / (1 - upsilon)) * aux1 * aux2;

% Apply collateral constraint
k_star = min(k_unc, lambda * a);

% Compute total labor input from the labor FOC
exp1   = 1 / (1 - (1 - gamma) * upsilon);
n_bar  = (theta * (1 - gamma) * upsilon / w)^exp1 * k_star^(gamma * upsilon * exp1);

% Check whether hired labor is strictly positive

if n_bar > lbar
    % Constraint n >= 0 is NOT binding
    k = k_star;
    n = n_bar - lbar;

else
    % Case 2: Re-solve imposing n = 0, so only own labor is used

    aux1  = (gamma / (r + delta))^(1 / (1 - gamma * upsilon));
    aux2  = lbar^((upsilon * (1 - gamma)) / (1 - gamma * upsilon));
    k_unc = (theta * upsilon)^(1 / (1 - gamma * upsilon)) * aux1 * aux2;

    % Apply collateral constraint
    k_star = min(k_unc, lambda * a);

    k = k_star;
    n = 0;
    n_bar = lbar;
end

%% Output (same formulas in both cases)

output = theta * (k^gamma * n_bar^(1 - gamma))^upsilon;

end %end function
