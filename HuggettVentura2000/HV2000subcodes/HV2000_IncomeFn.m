function income=HV2000_IncomeFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)
% Income as defined in HV2000 Section 5.1: labor earnings + capital income
% (before tax) + value of all transfers received (lump-sum T plus SS benefit).
% Used to compute Income Gini and top X% shares.

if agej<Jr
    earnings=z*ybar_j*w;
    benefit=0;
else
    earnings=0;
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

income=earnings+r*a+T+benefit;

end
