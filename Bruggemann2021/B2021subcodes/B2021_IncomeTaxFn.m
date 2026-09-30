function TaxI=B2021_IncomeTaxFn(I,d,tau_s, ybar, tau_i_r1,tau_i_r2,tau_i_r3,tau_i_r4,tau_i_r5,tau_i_r6, tau_i_t1, tau_i_t2, tau_i_t3, tau_i_t4, tau_i_t5, tau_i_t6)
% B2021_IncomeTaxFn gives progressive income tax schedule.
% Applies the deduction, bracketed federal income tax rates, and the
% proportional state and local tax tau_s to taxable income.

% Taxable income
Id = I-d*ybar; % Income minus deduction (scaled by ybar)
Id = max(Id,0); % This is done before we calculate the taxes
% this eqn max(Id,0) is not in paper, but is in B2021 codes, 
% function tax(yy) and lines such as
% vectaxe(counter)=tax(grossinc(counter))+taubal*(max(grossinc(counter)-deduc,0.0_DP)) 

TaxI = 0; % keep the gpu happy
% Calculate the tax based on the applicable rates and thresholds
% Note: following hardcodes that tau_i_t1=0
% Note that all the thresholds are scaled by ybar
if tau_i_t1*ybar < Id && Id <= tau_i_t2*ybar
    TaxI = Id * tau_i_r1;
elseif tau_i_t2*ybar < Id && Id <= tau_i_t3*ybar
    TaxI = tau_i_t2*ybar * tau_i_r1 + (Id - tau_i_t2*ybar) * tau_i_r2;
elseif  tau_i_t3*ybar < Id && Id <= tau_i_t4*ybar
    TaxI = tau_i_t2*ybar * tau_i_r1 + (tau_i_t3*ybar - tau_i_t2*ybar) * tau_i_r2 + (Id - tau_i_t3*ybar) * tau_i_r3;
elseif  tau_i_t4*ybar < Id && Id <= tau_i_t5*ybar
    TaxI = tau_i_t2*ybar * tau_i_r1 + (tau_i_t3*ybar - tau_i_t2*ybar) * tau_i_r2 + (tau_i_t4*ybar - tau_i_t3*ybar) * tau_i_r3 + (Id - tau_i_t4*ybar) * tau_i_r4;
elseif  tau_i_t5*ybar < Id && Id <= tau_i_t6*ybar
    TaxI = tau_i_t2*ybar * tau_i_r1 + (tau_i_t3*ybar - tau_i_t2*ybar) * tau_i_r2 + (tau_i_t4*ybar - tau_i_t3*ybar) * tau_i_r3 + (tau_i_t5*ybar - tau_i_t4*ybar) * tau_i_r4 + (Id - tau_i_t5*ybar) * tau_i_r5;
elseif tau_i_t6*ybar < Id
    TaxI = tau_i_t2*ybar * tau_i_r1 + (tau_i_t3*ybar - tau_i_t2*ybar) * tau_i_r2 + (tau_i_t4*ybar - tau_i_t3*ybar) * tau_i_r3 + (tau_i_t5*ybar - tau_i_t4*ybar) * tau_i_r4 + (tau_i_t6*ybar - tau_i_t5*ybar) * tau_i_r5 + (Id - tau_i_t6*ybar) * tau_i_r6;
end

TaxI = TaxI+tau_s*Id;

end %end function
