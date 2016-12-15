function F=HubbardSkinnerZeldes1994_ReturnFn(aprime,a,W_z1,M_z2,gamma,r,Cbar,DeterministicWj, DeterministicMj) %i

%i: education type
%jj: age (index 1:J; indicates ages 21 to 100)

Wj=DeterministicWj+exp(W_z1);
Mj=DeterministicMj+exp(M_z2);

c=(1+r)*a+Wj-Mj-aprime;

F=-Inf; %(l by aprime)

if c>Cbar
    F=(c.^(1-gamma)-1)/(1-gamma);
else
    F=(Cbar.^(1-gamma)-1)/(1-gamma);
end


end