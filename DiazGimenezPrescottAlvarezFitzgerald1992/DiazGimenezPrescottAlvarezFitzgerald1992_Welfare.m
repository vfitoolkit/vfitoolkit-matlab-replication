function [ValueFnofPublicandPrivateConsumption,v1,v2,W]=DiazGimenezPrescottAlvarezFitzgerald1992_Welfare(Params,V,n_A,n_K,n_s,n_z,n_d,n_a,n_sz,sigma_sz,w_sz,pi_sz,A_grid,K_grid,d_grid, a_grid,sz_grid,StationaryDist, Policy)
% Calculates V, v1, v2, and W from formulae on page 553


% h(s,z): lifetime human capital wealth
hdist=Inf;
h=zeros(n_s,n_z);
while hdist>10^(-9)
    hold=h;
    h=sigma_sz.*w_sz+Params.beta*sigma_sz.*(pi_sz*hold); % might need to be pi_sz'

    hdist=sum(abs(h-hold));    
end

% Total Wealth of Household
W=zeros(n_A,n_K,n_s,n_z);
for i1=1:n_A
    for i2=1:n_K
        for i3=1:n_s
            for i4=1:n_z
                W(i1,i2,i3,i4)=A_grid(i1)+K_grid(i2)+h(i3,i4);
            end
        end
    end
end

% v1 is the value function computed earlier

% v2: value of public consumption (note: this code only works for n_z=1)
SSvalueParamNames=struct();
SSvalueParamNames(1).Names={'phi'};
SSvalueParamNames(2).Names={'eta_d','eta_l'};
SSvalueParamNames(3).Names={'w1','w2'};
SSvalueParamNames(4).Names={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};
SSvalueParamNames(5).Names={};
SSvaluesFn_1 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi) a2prime_val*(a2prime_val-a2_val>0)-phi*a2_val*(a2prime_val-a2_val<0); % Investment (Housing investment)
SSvaluesFn_2 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,eta_d,eta_l) eta_d*a1_val*(a1_val>0)-eta_l*a1_val*(a1_val<0); % Banking Value Added
SSvaluesFn_3 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,w1,w2) w1*(s_val==1)*d_val+w2*(s_val==2)*d_val; % Wage income
SSvaluesFn_4 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) DiazGimenezPrescottAlvarezFitzgerald1992_ConsumptionFn(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment); % Goods Consumption of Homeowners
SSvaluesFn_5 = @(d_val,a1prime_val,a2prime_val,a1_val,a2_val,s_val,z_val) a2_val; % Capital (Housing Stock)
SSvaluesFn={SSvaluesFn_1, SSvaluesFn_2, SSvaluesFn_3, SSvaluesFn_4, SSvaluesFn_5};
SSvalues_AggVars3=SSvalues_AggVars_Case1(StationaryDist, Policy, SSvaluesFn, Params, SSvalueParamNames, n_d, n_a, n_sz, d_grid, a_grid,sz_grid, 2);
g=8*(SSvalues_AggVars3(3)-(SSvalues_AggVars3(4)+SSvalues_AggVars3(1)+SSvalues_AggVars3(2)+Params.mew*SSvalues_AggVars3(5)));

v2dist=Inf;
v2=zeros(n_s,n_z);
while v2dist>10^(-9)
    v2old=v2;
    v2=sigma_sz*Params.delta_g*(g^(Params.alpha*(1-Params.psi)))/(1-Params.psi)+Params.beta*sigma_sz.*(pi_sz*v2old);
    
    v2dist=sum(abs(v2-v2old));
end
v2=real(v2); % for some reason matlab is making them complex numbers even though imaginary component is 0
v2=gather(v2);

% V= v1+v2
v1=gather(V);
ValueFnofPublicandPrivateConsumption=zeros(n_A,n_K,n_s,n_z);
for i1=1:n_A
    for i2=1:n_K
        for i3=1:n_s
            for i4=1:n_z
                ValueFnofPublicandPrivateConsumption(i1,i2,i3,i4)=v1(i1,i2,i3,i4)+v2(i3,i4);
            end
        end
    end
end

end