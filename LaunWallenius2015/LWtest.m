% LWtest.m -- standalone test discriminating Reading A vs Reading B for the
% LW2015 health-transition Markov chain.
%
% In LW2015 health is purely exogenous (no mortality, no behavioral
% feedback), so the unconditional health distribution at any age is just
% the initial distribution propagated through pi_z_J. LW2015 p.132 reports
% the model's distribution at age 70 as (42, 25, 13, 9, 11) % across
% (very good, good, fair, bad, very bad). This file builds pi_z_J under
% both readings, propagates the initial distribution P(h=5)=1 at j=1 to
% j=46 (real age 70), and prints both side-by-side with the paper target.
% Whichever reading reproduces the paper's numbers is the one LW2015 used.
%
% Reading A (independent Bernoulli draws): each shock fires independently
%   each period; the joint distribution has 8 outcomes; net change in h is
%   the algebraic sum capped to [1,5].
% Reading B (mutually exclusive with renormalization): the three shock
%   probabilities are treated as competing event probabilities; when their
%   sum exceeds 1 they are rescaled to sum to 1 (no "stay" mass) and
%   p_pos is also forced to 0 at h=5 (since it's a no-op at the max).

clear

agejshifter=24;
N_j=80-agejshifter;
p_pos_young=[0.6 0.7 0.8 0.9 1.0];

%% Build pi_z_J under Reading A (independent Bernoulli draws)
pi_z_J_A=zeros(5,5,N_j);
for jj=1:N_j
    age_frac=(jj-1)/(N_j-1);
    for hh=1:5
        if hh==1
            p_l=0.6; p_h=0.1;
        else
            p_l=0.1+0.5*age_frac;
            p_h=0.01;
        end
        p_i=p_pos_young(hh)-(p_pos_young(hh)-0.5)*age_frac;

        for s_l=0:1
            for s_h=0:1
                for s_i=0:1
                    prob=(s_l*p_l+(1-s_l)*(1-p_l)) ...
                        *(s_h*p_h+(1-s_h)*(1-p_h)) ...
                        *(s_i*p_i+(1-s_i)*(1-p_i));
                    delta=-s_l-3*s_h+s_i;
                    h_prime=max(1,min(5,hh+delta));
                    pi_z_J_A(hh,h_prime,jj)=pi_z_J_A(hh,h_prime,jj)+prob;
                end
            end
        end
    end
end

%% Build pi_z_J under Reading B (mutually exclusive + renormalize on overshoot)
pi_z_J_B=zeros(5,5,N_j);
for jj=1:N_j
    age_frac=(jj-1)/(N_j-1);
    p_l_normal=0.1+0.5*age_frac;
    p_h_normal=0.01;
    for hh=1:5
        if hh==1
            p_neg_lo=0;
            p_neg_hi=0;
            p_pos=p_pos_young(1)-(p_pos_young(1)-0.5)*age_frac;
        elseif hh==5
            p_neg_lo=p_l_normal;
            p_neg_hi=p_h_normal;
            p_pos=0;
        else
            p_neg_lo=p_l_normal;
            p_neg_hi=p_h_normal;
            p_pos=p_pos_young(hh)-(p_pos_young(hh)-0.5)*age_frac;
        end
        total=p_neg_lo+p_neg_hi+p_pos;
        if total>1
            p_neg_lo=p_neg_lo/total;
            p_neg_hi=p_neg_hi/total;
            p_pos=p_pos/total;
            p_stay=0;
        else
            p_stay=1-total;
        end
        pi_z_J_B(hh,max(1,hh-1),jj)=pi_z_J_B(hh,max(1,hh-1),jj)+p_neg_lo;
        pi_z_J_B(hh,max(1,hh-3),jj)=pi_z_J_B(hh,max(1,hh-3),jj)+p_neg_hi;
        pi_z_J_B(hh,min(5,hh+1),jj)=pi_z_J_B(hh,min(5,hh+1),jj)+p_pos;
        pi_z_J_B(hh,hh,jj)=pi_z_J_B(hh,hh,jj)+p_stay;
    end
end

%% Propagate initial dist P(h=5)=1 from j=1 to age 70 (j=46) under each
% z_grid order is h=1..5 = (very bad, bad, fair, good, very good)
h_dist_A=[0 0 0 0 1];
h_dist_B=[0 0 0 0 1];
for jj=1:45
    h_dist_A=h_dist_A*pi_z_J_A(:,:,jj);
    h_dist_B=h_dist_B*pi_z_J_B(:,:,jj);
end

%% Report -- order matches paper p.132 (very good, good, fair, bad, very bad)
A_pct=100*[h_dist_A(5) h_dist_A(4) h_dist_A(3) h_dist_A(2) h_dist_A(1)];
B_pct=100*[h_dist_B(5) h_dist_B(4) h_dist_B(3) h_dist_B(2) h_dist_B(1)];
target=[42 25 13 9 11];

fprintf('\n');
fprintf('Health distribution at real age 70 (model age j=46), %%\n');
fprintf('Order: very good | good | fair | bad | very bad\n\n');
fprintf('  LW2015 p.132 reported model:  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n',target);
fprintf('  Reading A (indep draws):      %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n',A_pct);
fprintf('  Reading B (mut excl + renorm):%5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n\n',B_pct);

% Quick L1 distance summary
fprintf('  L1 distance to paper:  Reading A = %.1f,  Reading B = %.1f\n', ...
    sum(abs(A_pct-target)), sum(abs(B_pct-target)));
fprintf('\n');
