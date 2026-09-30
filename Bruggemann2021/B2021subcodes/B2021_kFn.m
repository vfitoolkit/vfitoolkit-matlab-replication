function k=B2021_kFn(l,e,aprime,a,age,eta,theta, r,w,lambda,delta,gamma,upsilon,lbar)
% B2021_kFn Entrepreneurial capital demand.
% Returns zero for non-entrepreneurs and retirees, and otherwise recovers
% entrepreneurial capital from f_solve_entre with own labor fixed at lbar.

k=0;

if e==1 && age==1 % Entrepreneurs' problem
    [k,~,~]=B2021_EntrepreneurStaticProblem(a,theta,r,w,lambda,delta,gamma,upsilon,lbar);
end

end