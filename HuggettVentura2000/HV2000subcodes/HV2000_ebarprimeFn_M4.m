function ebarprime=HV2000_ebarprimeFn_M4(d,ebar,z,e,agej,Jr,w,ybar_j,ybar_mean,g)
% Model 4 ebar law of motion (experienceassetze signature aprime(d,a,z,e)).
% Earnings include the iid e shock: min(z*e*ybar_j*w, ebar_max).

ebarprime=0;

if agej<Jr
    ebar_max=2.47*w*ybar_mean;
    earnings=min(z*e*ybar_j*w,ebar_max);
    ebarprime=(ebar*(agej-1)+earnings)/agej;
else
    ebarprime=ebar/(1+g);
end

end
