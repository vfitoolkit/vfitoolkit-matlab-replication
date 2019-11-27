function [e,pi_s,gammastar,gammastarfull]=CastanedaDiazGimenezRiosRull2003_Create_Exog_Shock(Params)
% Creates the transition matrix on the exogenous shocks.

% Note, there are both gammastar and gammastarfull

e=gpuArray([1,Params.e2,Params.e3,Params.e4, zeros(1,Params.J)]); %e(s)

Gamma_ee_11=1-Params.Gamma_ee_12-Params.Gamma_ee_13-Params.Gamma_ee_14-Params.p_eg;
Gamma_ee_22=1-Params.Gamma_ee_21-Params.Gamma_ee_23-Params.Gamma_ee_24-Params.p_eg;
Gamma_ee_33=1-Params.Gamma_ee_31-Params.Gamma_ee_32-Params.Gamma_ee_34-Params.p_eg;
Gamma_ee_44=1-Params.Gamma_ee_41-Params.Gamma_ee_42-Params.Gamma_ee_43-Params.p_eg;

Gamma_ee=gpuArray([Gamma_ee_11, Params.Gamma_ee_12, Params.Gamma_ee_13, Params.Gamma_ee_14; Params.Gamma_ee_21, Gamma_ee_22, Params.Gamma_ee_23, Params.Gamma_ee_24; Params.Gamma_ee_31, Params.Gamma_ee_32, Gamma_ee_33, Params.Gamma_ee_34; Params.Gamma_ee_41, Params.Gamma_ee_42, Params.Gamma_ee_43, Gamma_ee_44]);

Gamma_ee=Gamma_ee./(sum(Gamma_ee,2)*ones(1,Params.J,'gpuArray')); %This is a normalization of Gamma_ee into a probability matrix
% Gamma_ee=Gamma_ee./(1-p_eg); %this would be equivalent to the normalization in the line above

% Calculate the gammastar's
gammastar=ones(1,Params.J)/Params.J;
for ii=1:5000
    gammastar=gammastar*Gamma_ee;
end

%Now we calculate all of the points for Gamma_re
%Step 1
p_51_Step1=gammastar(1)+Params.phi1*gammastar(2)+Params.phi1^2*gammastar(3)+Params.phi1^3*gammastar(4);
p_52_Step1=(1-Params.phi1)*(gammastar(2)+Params.phi1*gammastar(3)+Params.phi1^2*gammastar(4));
p_53_Step1=(1-Params.phi1)*(gammastar(3)+Params.phi1*gammastar(4));
p_54_Step1=(1-Params.phi1)*gammastar(4);
p_61_Step1=(1-Params.phi1)*gammastar(1);
p_62_Step1=Params.phi1*gammastar(1)+gammastar(2)+Params.phi1*gammastar(3)+Params.phi1^2*gammastar(4);
p_63_Step1=(1-Params.phi1)*(gammastar(3)+Params.phi1*gammastar(4));
p_64_Step1=(1-Params.phi1)*gammastar(4);
p_71_Step1=(1-Params.phi1)*gammastar(1);
p_72_Step1=(1-Params.phi1)*(Params.phi1*gammastar(1)+gammastar(2));
p_73_Step1=Params.phi1^2*gammastar(1)+Params.phi1*gammastar(2)+gammastar(3)+Params.phi1*gammastar(4);
p_74_Step1=(1-Params.phi1)*gammastar(4);
p_81_Step1=(1-Params.phi1)*gammastar(1);
p_82_Step1=(1-Params.phi1)*(Params.phi1*gammastar(1)+gammastar(2));
p_83_Step1=(1-Params.phi1)*(Params.phi1^2*gammastar(1)+Params.phi1*gammastar(2)+gammastar(3));
p_84_Step1=Params.phi1^3*gammastar(1)+Params.phi1^2*gammastar(2)+Params.phi1*gammastar(3)+gammastar(4);
%Step 2
p_51=p_51_Step1+Params.phi2*p_52_Step1+Params.phi2^2*p_53_Step1+Params.phi2^3*p_54_Step1;
p_61=p_61_Step1+Params.phi2*p_62_Step1+Params.phi2^2*p_63_Step1+Params.phi2^3*p_64_Step1;
p_71=p_71_Step1+Params.phi2*p_72_Step1+Params.phi2^2*p_73_Step1+Params.phi2^3*p_74_Step1;
p_81=p_81_Step1+Params.phi2*p_82_Step1+Params.phi2^2*p_83_Step1+Params.phi2^3*p_84_Step1;
p_52=(1-Params.phi2)*(p_52_Step1+Params.phi2*p_53_Step1+Params.phi2^2*p_54_Step1);
p_62=(1-Params.phi2)*(p_62_Step1+Params.phi2*p_63_Step1+Params.phi2^2*p_64_Step1);
p_72=(1-Params.phi2)*(p_72_Step1+Params.phi2*p_73_Step1+Params.phi2^2*p_74_Step1);
p_82=(1-Params.phi2)*(p_82_Step1+Params.phi2*p_83_Step1+Params.phi2^2*p_84_Step1);
p_53=(1-Params.phi2)*(p_53_Step1+Params.phi2*p_54_Step1);
p_63=(1-Params.phi2)*(p_63_Step1+Params.phi2*p_64_Step1);
p_73=(1-Params.phi2)*(p_73_Step1+Params.phi2*p_74_Step1);
p_83=(1-Params.phi2)*(p_83_Step1+Params.phi2*p_84_Step1);
p_54=(1-Params.phi2)*p_54_Step1;
p_64=(1-Params.phi2)*p_64_Step1;
p_74=(1-Params.phi2)*p_74_Step1;
p_84=(1-Params.phi2)*p_84_Step1;
%Now put these into the matrix
Gamma_re=gpuArray([p_51, p_52, p_53, p_54; p_61, p_62, p_63, p_64; p_71, p_72, p_73, p_74; p_81, p_82, p_83, p_84]);

Gamma_re=Gamma_re./(sum(Gamma_re,2)*ones(1,Params.J,'gpuArray')); %This is a normalization of Gamma_re into a probability matrix

pi_s=[Gamma_ee.*(1-Params.p_eg), diag(Params.p_eg*ones(Params.J,1,'gpuArray')); Gamma_re.*(1-Params.p_gg), diag(Params.p_gg*ones(Params.J,1,'gpuArray'))];  %transmatix is (z,zprime) %dim N_s-by-N_s (s by sprime)
% Make doubly sure that rows sum to exactly one (had a rounding error on
% very rare occasions, which was only apparent when solving thousands of
% times as part of calibration codes.
if min(min(pi_s))>=0
    fix_pi_z_counter=0;
    while (min(sum(pi_s,2))~=1 || max(sum(pi_s,2))~=1) && fix_pi_z_counter<10
        if min(abs(sum(pi_s,2)-1)) < 10^(-12)
            fprintf('adjusting pi_s so rows sum to 1 (first attempt) \n')
            [min(sum(pi_s,2)),max(sum(pi_s,2))]
            [min(sum(pi_s,2))-1,max(sum(pi_s,2))-1]
            pi_s=pi_s./sum(pi_s,2);
        end
        % If that does not fix it then try adjusting the diagonals
        if min(abs(sum(pi_s,2)-1)) < 10^(-12)
            fprintf('adjusting pi_s so rows sum to 1 (second attempt) \n')
            [min(sum(pi_s,2)),max(sum(pi_s,2))]
            [min(sum(pi_s,2))-1,max(sum(pi_s,2))-1]
            pi_s=pi_s-(eye(size(pi_s)).*(sum(pi_s,2)-1));
        end
        fix_pi_z_counter=fix_pi_z_counter+1;
%         if min(abs(sum(pi_s,2)-1)) < 10^(-12)
%             fprintf('adjusting pi_s so rows sum to 1 (third attempt) \n')
%             [min(sum(pi_s,2)),max(sum(pi_s,2))]
%             [min(sum(pi_s,2))-1,max(sum(pi_s,2))-1]
%             pi_s=pi_s./sum(pi_s,2);
%         end
    end
end

gammastarfull=ones(1,2*Params.J)/(2*Params.J);
for ii=1:5000
    gammastarfull=gammastarfull*pi_s;
end


end