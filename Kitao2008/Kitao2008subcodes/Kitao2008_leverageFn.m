function leverage=Kitao2008_leverageFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)
% Leverage ratio (k-a)/k of the entrepreneur, that is borrowing as a fraction of the investment k.
% Returns 0 if worker, or if the entrepreneur does not borrow.

leverage=0;
k=0; % just to make GPU happy

if e==1 % Entrepreneurs' problem (for workers, k stays at zero)
    % The entrepreneurs production problem is static
    [k,n_unused,y_unused]=Kitao2008_StaticEntrepreneurProblem(a,eta,r,w,d,phi,delta,theta,upsilon1,upsilon2,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);
end

% Borrowing is k-a, and the leverage ratio is that as a fraction of the investment k
if k>a
    leverage=(k-a)/k;
end
% Otherwise leverage is zero (the entrepreneur is not a borrower)
%
% Note on the definition: Kitao (2008) does not say how the "avg. lev. ratio (%)" column of her
% Table 5 is built, and this replication originally used (k-a)/a, borrowing as a fraction of own
% assets. That cannot be what she reports. An entrepreneur at the collateral constraint has
% k=(1+d)a, so (k-a)/a=d=0.5 exactly; since essentially every theta4 entrepreneur is constrained,
% the column came out at 49.7% against the 32.6% Kitao reports, and no amount of widening the asset
% grid moved it. Under (k-a)/k the same constrained entrepreneur gives 0.5a/1.5a=1/3 instead.
% A group with a constrained share p averages 0.5p under the old definition and p/3 under the new,
% so the whole column can be mapped across:
%    theta2   3.0% -> 2.0%    (Kitao 2.4%)
%    theta3  22.4% -> 14.9%   (Kitao 15.2%)
%    theta4  49.7% -> 33.1%   (Kitao 32.6%)
%    all     32.4% -> 21.6%   (Kitao 21.1%)
% All four rows line up, so (k-a)/k is taken to be the definition Kitao (2008) uses.

end
