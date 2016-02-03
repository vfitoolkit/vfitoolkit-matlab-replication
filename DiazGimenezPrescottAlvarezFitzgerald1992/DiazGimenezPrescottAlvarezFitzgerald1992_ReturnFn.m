function F=DiazGimenezPrescottAlvarezFitzgerald1992_ReturnFn(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment)
% A is called a by Diaz-Gimenez, Prescott, Alvarez & Fitzgerald (1992) [renamed so can use VFI toolkit standard notation of 'a' being endogenous state vector]

F=-Inf;
C=nan;

%firstly calculate the values of x_d, x_s
x_d=0;
x_s=0;
if K_val>Kprime_val
%     x_d=0;
    x_s=phi*(K_val-Kprime_val);
elseif K_val<Kprime_val
    x_d=Kprime_val-K_val;
%     x_s=0;
end


%secondly calculate the values of d & l
%if Aprime_val==0;
d=0;
l=0;
if Aprime_val>0
    d=e*Aprime_val;
%     l=0;
elseif Aprime_val<0
%     d=0;
    l=-Aprime_val*e; %the minus is because Aprime_val itself is negative
end

pension=0;
if A_val==0 && K_val==1 && s_val==3 %retired, no assets, no house
    pension=omega;
end

w=nan;
if s_val==1
    w=w1;
elseif s_val==2
    w=w2;
elseif s_val==3
    w=w3;
elseif s_val==4
    w=w4;
end
if Experiment==2
    if z_val==1
        w=w*1.02;
    elseif z_val
        w=w*0.98;
    end
end

C=A_val+(1-theta)*(w*N_val+d*i_d-l*i_l) + x_s-x_d +l-d+pension-mew*Kprime_val;

if C>0 %If this is not met, then just leave Fmatrix as -Inf
    if s_val==1 || s_val==2
        if K_val==0
            F=( ( C^alpha * (1-alpha_k/alpha)^(alpha-alpha_k) * (alpha_k/(gamma*alpha))^alpha_k * (tau-N_val)^(1-alpha) )^(1-psi) )/(1-psi);
        end
        if K_val>0 %ie. if K=3
            F=( ( C^(alpha-alpha_k) * K_val^alpha_k * (tau-N_val)^(1-alpha) )^(1-psi) )/(1-psi);
        end
    end
    if s_val==3
        F=delta_r* ((C^alpha)^(1-psi)) /(1-psi);
    end
    %notice that is s==4, then return will be zero, so no need to calculate
    %anything as we already set this when Fmatrix is first created
end

end