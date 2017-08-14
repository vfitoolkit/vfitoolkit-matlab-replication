function c=HubbardSkinnerZeldes1994_ConsumptionFn(aprime,a,W_z1,M_z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj) %i

%i: education type
%jj: age (index 1:J; indicates ages 21 to 100)


% Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% "This specification ensures that when we compare the certainty case with
% the earnings uncertainty case, we hold the age-conditional means of earnings constant."
Wj=exp(log(DeterministicWj)-0.5*w_sigmasqu+W_z1);
Mj=exp(DeterministicMj+M_z2);

TR=max(Cbar+Mj-(1+r)*a-Wj,0);

c=(1+r)*a+Wj-Mj+TR-aprime;

end