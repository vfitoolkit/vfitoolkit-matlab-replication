function n=B2021_nFn(l,e,aprime,a,age,eta,theta, r,w,lambda,delta,gamma,upsilon,lbar)
% B2021_nFn Entrepreneurial hired labor demand.
% Returns hired labor only, excluding the entrepreneur's own labor lbar,
% and is zero for workers and retirees.

n=0; 

if e==1 && age==1 % Entrepreneurs' problem

    [~,n,~]=B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar);

end 

end %end function
