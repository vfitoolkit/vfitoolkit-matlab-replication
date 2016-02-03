function c=Imrohoroglu1989_ConsFn(aprime_val, a_val, s_val,z_val,r_l,r_b,y,theta,sigma)

% s_val=rem(sz_val-1,2)+1;

%F=-Inf; %model: (aprime,a,s)
c=0;

if aprime_val>=0 %If lending
    r=r_l;
else %If borrowing
    r=r_b;
end

earnings=0;
if s_val==1 %If employed
    earnings=y;
elseif s_val==2 %If unemployed
    earnings=theta*y;
end

c=a_val-aprime_val/(1+r)+earnings;

% if c>0
%     F=(c^(1-sigma))/(1-sigma);
% end

end