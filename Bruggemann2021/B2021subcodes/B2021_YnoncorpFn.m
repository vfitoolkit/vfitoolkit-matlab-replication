function output=B2021_YnoncorpFn(l,e,aprime,a,age,eta,theta,r,w,lambda,delta,gamma,upsilon,lbar)
% B2021_YnoncorpFn Entrepreneurial sector output.
% Returns non-corporate output for young entrepreneurs and zero otherwise,
% with firm choices recovered from f_solve_entre using lbar.

output=0;

if e==1 && age==1 % Entrepreneurs' problem

    [~,~,output]=B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar);
    
end %end function












end
