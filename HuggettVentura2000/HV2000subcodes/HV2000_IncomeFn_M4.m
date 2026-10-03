function income=HV2000_IncomeFn_M4(d,aprime,a,ebar,z,e,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)
% Model 4 income (mirrors HV2000_IncomeFn; earnings include iid e shock).

if agej<Jr
    earnings=z*e*ybar_j*w;
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
