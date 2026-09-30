function I=B2021_IncomeFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar,pension)
% B2021_IncomeFn gives household pre-tax income.
% Returns taxable income for workers, entrepreneurs, and retirees, with
% entrepreneurial income recovered from B2021_EntrepreneurStaticProblem.
%
% Note: omits any lumpsum that might be received

I=0; % B2021 denotes this y, but I am using I in line with the contents of the ReturnFn

if e==0 && age==1 % Young workers' problem
    
    I=w*l*eta+r*a; % income, y^w

elseif e==1 && age==1 % Entrepreneurs' problem

    [k,n,output]=B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar);
    
    I=output-delta*k-r*(k-a)-w*n; % Income, y^e
            
elseif age==2 % Retiree's problem (for both e=0 & e=1)

    I=r*a+pension; % income, y^r

end

end %end function
