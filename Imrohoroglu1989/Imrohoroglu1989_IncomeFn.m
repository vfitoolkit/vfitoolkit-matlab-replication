function income=Imrohoroglu1989_IncomeFn(aprime_val, a_val, s_val,z_val,r_l,r_b,y,theta,sigma)

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

income=a_val*r+earnings;

end