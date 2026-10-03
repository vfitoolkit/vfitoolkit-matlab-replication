function c=HV2000_ConsumptionFn_M4(d,aprime,a,ebar,z,e,agej,Jr,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,haveSS)
% Model 4 consumption (mirrors HV2000_ConsumptionFn but with iid e shock).

if agej<Jr
    earnings=z*e*ybar_j*w;
else
    earnings=0;
end

if agej<Jr
    benefit=0;
else
    b_bend1=0.20*w*ybar_mean;
    b_bend2=1.24*w*ybar_mean;
    if ebar<=b_bend1
        b_er=0.90*ebar;
    elseif ebar<=b_bend2
        b_er=0.90*b_bend1+0.32*(ebar-b_bend1);
    else
        b_er=0.90*b_bend1+0.32*(b_bend2-b_bend1)+0.15*(ebar-b_bend2);
    end
    benefit=haveSS*(b_common+b_er/((1+g)^(agej-Jr)));
end

c=a*(1+r*(1-tau))+(1-tau-theta)*earnings+T+benefit-(1+g)*aprime;

end
