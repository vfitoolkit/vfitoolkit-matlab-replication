function [k,n,output]=Kitao2008_StaticEntrepreneurProblem(a,eta,r,w,d,phi,delta,theta,upsilon1,upsilon2,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)
% The entrepreneurs production problem is static.
% See http://discourse.vfitoolkit.com/t/on-solving-models-with-entrepreneurs/296
%
% Solves  max_{k,n} theta*k^upsilon1*n^upsilon2 - w*max(n-eta,0) - (rbar+delta)*k   s.t. k<=(1+d)*a
%
% Two complications relative to the textbook decreasing-returns-to-scale problem:
% (i)  there is a different r for borrowing and for saving, rbar=r+phi if k>a and rbar=r if k<=a
% (ii) the entrepreneur uses their own labor endowment eta in their own project at zero cost, and
%      only buys labor above eta at the wage w (Kitao (2008), eqn (7) and footnote 6). So the
%      marginal cost of labor is zero for n<eta and w for n>eta. Hence the optimal n is never
%      below eta (except when it is optimal to produce nothing at all), and when the n=eta kink
%      binds the first-order condition for k changes (labor is then fixed at eta, not chosen).
% (iii) ONLY UNDER taxincome=2 the capital choice has to be made after tax, not before. Eqn (7) of
%      Kitao (2008) has the tax liability T(I) inside the max, so k and n maximise profit net of
%      tax. The objective is f+(1-delta)k-(1+rbar)(k-a)-w*max(n-eta,0)-T, whose non-tax part is
%      just I+a with a fixed, so the entrepreneur maximises I-T. Whether that has the same argmax
%      as maximising I alone depends on the tax regime:
%        taxincome=1: T=T(I), and I-T(I) is increasing in I whenever T'<1, so the argmax is
%                     unchanged and maximising pre-tax income is correct.
%        taxincome=3: I=IE1+IE2 with IE2=r*a independent of k and n, so I-T=(1-tau_E1)*IE1 plus
%                     terms that do not involve k or n. Again the argmax is unchanged.
%        taxincome=2: I-T=(1-tau_k)*IK + Itilde-G(Itilde)-tau_I*Itilde, and BOTH IK=r*(a-k) and
%                     Itilde move with k, at different marginal rates. The argmax genuinely
%                     differs, and this is the margin Kitao (2008) describes in Section 5.1:
%                     "entrepreneurs find saving a more attractive use of assets relative to
%                     entrepreneurial investment" when the capital income tax is low.
%      For k<a the r terms cancel out of Itilde, leaving Itilde=f(k,n,theta)-delta*k-w*max(n-eta,0),
%      so dItilde/dk=f_k-delta and the first-order condition is
%           [1-mtr(Itilde)]*(f_k-delta) = (1-tau_k)*r,   i.e.   f_k = delta + r_eff
%      with r_eff = r*(1-tau_k)/(1-mtr) and mtr = G'(Itilde)+tau_I the marginal rate of the
%      progressive schedule. So the whole change is that the entrepreneur discounts at an effective
%      cost of capital r_eff rather than r: a low tau_k raises r_eff and crowds capital out of the
%      business, a high tau_k lowers it and draws capital in.
%      For k>a there is no riskless saving, IK=0, and the monotone argument applies again, so the
%      borrowing branch below (k_unc, based on r+phi) is already the after-tax solution in every
%      regime and is left alone.
%      mtr depends on Itilde which depends on k, so r_eff is implicit. It is solved by a fixed
%      three-step iteration rather than a convergence loop, so that the trip count is known at
%      compile time and the function stays gpuArray/arrayfun compilable.
onegg=1-upsilon1-upsilon2; % just to simplify the below forumlaes

k=0; % just to make GPU happy
n=0;
output=0;

if theta>0 && a>0 % otherwise the project produces nothing, and it is optimal to shut down (k=n=y=0)

    % Find capital used by entrepreneur if they borrow (so based on r+phi)
    k_unc=(upsilon1/(r+phi+delta))^((1-upsilon2)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
    if ((upsilon2*theta*(k_unc^upsilon1))/w)^(1/(1-upsilon2))<eta
        % Own labor endowment is not exhausted at the interior solution, so n=eta and the marginal
        % cost of labor is zero. FOC for k becomes upsilon1*theta*eta^upsilon2*k^(upsilon1-1)=rbar+delta
        k_unc=(upsilon1*theta*(eta^upsilon2)/(r+phi+delta))^(1/(1-upsilon1));
    end
    k=Inf; % establish k immediately before the if/else that sets it, to keep the GPU happy. Inf rather
    % than 0 so that if it ever fails to be overwritten it shows up loudly rather than looking like a
    % legitimate shut-down (k=0 is what a genuine shut-down returns, see the initialisation above)
    if k_unc>(1+d)*a % Collateral constraint
        k=(1+d)*a;
    else
        k=k_unc;
    end
    % If k<a, switch to using the r for savings
    if k<a
        k_r=Inf; % establish before the branches that set it, as for k above
        if taxincome==2
            % After-tax capital choice, see (iii) at the top of this file. Only this regime needs it.
            % Fixed three-step iteration on r_eff; each step is the same closed form as the else
            % branch below, just with r_eff in place of r.
            r_eff=r; % first pass uses the pre-tax cost of capital
            n_r=0; Itilde=0; mtr=0; % establish everything the loop assigns, to keep the GPU happy
            for reffiter=1:3
                k_r=(upsilon1/(r_eff+delta))^((1-upsilon2)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
                if ((upsilon2*theta*(k_r^upsilon1))/w)^(1/(1-upsilon2))<eta
                    k_r=(upsilon1*theta*(eta^upsilon2)/(r_eff+delta))^(1/(1-upsilon1));
                end
                % Itilde at this k. The r terms cancel for k<a, so this is just output net of
                % depreciation and the wage bill.
                n_r=((upsilon2*theta*(k_r^upsilon1))/w)^(1/(1-upsilon2));
                if n_r<eta
                    n_r=eta;
                end
                Itilde=theta*(k_r^upsilon1)*(n_r^upsilon2)-delta*k_r-w*(n_r-eta)*(n_r>eta);
                % Marginal rate of the progressive schedule plus the proportional rate. G(I) is
                % tau_a0*(I-(I^(-tau_a1)+tau_a2)^(-1/tau_a1)), so G'(I) is the expression below.
                mtr=tau_I;
                if Itilde>0
                    mtr=tau_a0*(1-(tau_a2+Itilde^(-tau_a1))^(-1/tau_a1-1)*Itilde^(-tau_a1-1))+tau_I;
                end
                r_eff=r*(1-tau_k)/(1-mtr);
            end
        else
            % taxincome=1 and taxincome=3: maximising pre-tax income has the same argmax, see (iii).
            % Same two branches as above, but based on r rather than r+phi
            k_r=(upsilon1/(r+delta))^((1-upsilon2)/onegg) * (upsilon2/w)^(upsilon2/onegg) * theta^(1/onegg);
            if ((upsilon2*theta*(k_r^upsilon1))/w)^(1/(1-upsilon2))<eta
                k_r=(upsilon1*theta*(eta^upsilon2)/(r+delta))^(1/(1-upsilon1));
            end
        end
        if k_r>a
            k=a;
        else
            k=k_r;
        end
    end

    % Now find labor (the first eta units are the entrepreneur's own, and so are free)
    n=Inf; % as for k above: establish n immediately before the assignments that set it
    n=((upsilon2*theta*(k^upsilon1))/w)^(1/(1-upsilon2));
    if n<eta
        n=eta;
    end

    % Now that we have the production decisions for k and n, rest is straightforward
    output=Inf; % as for k and n above: establish output immediately before the assignment that sets it
    output=theta*(k^upsilon1)*(n^upsilon2);

end

end
