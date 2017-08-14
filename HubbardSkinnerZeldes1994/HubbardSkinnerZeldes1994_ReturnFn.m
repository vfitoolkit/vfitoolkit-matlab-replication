function F=HubbardSkinnerZeldes1994_ReturnFn(aprime,a,W_z1,M_z2,gamma,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj) %i

%i: education type
%jj: age (index 1:J; indicates ages 21 to 100)


% Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% "This specification ensures that when we compare the certainty case with
% the earnings uncertainty case, we hold the age-conditional means of earnings constant."
Wj=exp(log(DeterministicWj)-0.5*w_sigmasqu+W_z1);
Mj=exp(DeterministicMj+M_z2);

TR=max(Cbar+Mj-(1+r)*a-Wj,0);

c=(1+r)*a+Wj-Mj+TR-aprime;

F=-Inf; %(l by aprime)

if c>0
    if gamma~=1
        F=(c.^(1-gamma)-1)/(1-gamma);
    else
        F=log(c);
    end
end

% Paper is unclear about whether people at the consumption floor Cbar
% actually have to consume Cbar, or if they just receive enough transfers
% to reach Cbar but then decide for themselves how much of that to consume
% and how much to save. I have interpreted the floor as the later, since
% this appears to be how the transfers function implies that it works.

end