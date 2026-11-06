% DFJ2010_FirstStageParams: interest rate, taxes, health transitions, survival, income, medical-expense process
% (DFJ2010's first-stage estimates, from their replication materials), for all ten permanent types.
% Script shared by DeNardiFrenchJones2010.m (exogenous medical expenses) and DFJ2010_EndogMedex.m; needs n_z, N_j,
% Names_i in the workspace. [Moved out of DeNardiFrenchJones2010.m on 2026-09-29, line for line.]

% Interest rate
Params.r=1/0.98-1; % mu_r=1/0.98-1 in singlsav15_wave8.gau line 222 ("Mean net ROR"); approx 0.0204

% Goverment transfers: cfloor, set with the other preference parameters above

% Income taxes (from lines 189-194 of their wealth20_Singles.c)
% Based on tax brackets and marginal tax rates for each bracket
Params.taxbracket1=6250;
Params.taxbracket2=40200;
Params.taxbracket3=68400;
Params.taxbracket4=93950;
Params.taxbracket5=148250;
Params.taxbracket6=284700;
Params.margtaxrate0=0.0765; % tax rate 0 is for < bracket 1
Params.margtaxrate1=0.2616; % tax rate 1 is for brackets 1-to-2
Params.margtaxrate2=0.4119;
Params.margtaxrate3=0.3499;
Params.margtaxrate4=0.3834;
Params.margtaxrate5=0.4360;
Params.margtaxrate6=0.4761; 

% Estate taxes
Params.tau_e=0.15; % estate tax rate (DFJ2010 call this tautilde); tauBeq=0.15 in singlsav15_wave8.gau line 220
Params.estateexemption=600000; % estate tax exemption (DFJ2010 call this xtilde); exBeq=6e5 in singlsav15_wave8.gau line 221
% These are not in paper, got them from DFJ2010 codes.
% wealth20_Singles.c line 1782-3:
% tauBeq = asstvecPtr[1]; /* Estate tax rate */
% exBeq  = asstvecPtr[2]; /* Estate tax exemption level */
% which is read from a file, line 394
% asstvecPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"asstvec.fmt"));
% "asstvec.fmt" itself is not in the replication materials (it is generated at
% runtime), but opt/gmmcrit25.g line 35 builds it as asstvec=cfloor|tauBeq|exBeq|mu_r|sigma_r
% and singlsav15_wave8.gau lines 217-228 set those values.

% Medical expenses 1/2. (from Table 1):
% Two stochastic components: zeta which is AR(1) markov, and xi which is iid
% Note: these later get 'scaled up' by the sigma_coeff
Params.rho_m=0.922; % autocorrelation, persistent part
Params.sigma_mepsilon=sqrt(0.050); % std dev of innovations, persistent part
Params.sigma_mxi=sqrt(0.665); % std dev of innovations, transitory part
% Footnote 3: We assume that medical expenses are log-normally distributed, 
% so the predicted level of medical expenses is exp (m + 1/2*sigma^2), where 
% m denotes predicted log medical expenses and sigma^2 denotes the variance of 
% the idiosyncratic shock psi. The ratio of the level of medical expenses 
% two standard deviations above the mean to average medical expenses is
%    exp(m+2*sigma)/exp(m+1/2*sigma^2) =exp(2*sigma-1/2*sigma^2)
% which =6.8 if sigma=sqrt(2.53)
% Added Note: I'm sure why they use 2.53 as being a relevant number?

% Medical expenses 2/2.
% Two deterministic coefficients, m(g,h,I,j) and sigma(g,h,I,j) from eqn (6) of DFJ2010
% Call them m_coeff and sigma_coeff, and due to the dependence on
% health status I create four parameters (m_coeff_healthgood, m_coeff_healthbad, 
% sigma_coeff_healthgood, sigma_coeff_healthbad) each of which is
% age-and-ptype dependent. Then in ReturnFn use health status in
% if-statement to decide which of _healthgood and _healthbad to use).
% sigma_coeff scales the standard deviation of the stochastic component
% This is setup below.

% Medical expenses in the model are exp(m+sigma*psi), psi=zeta+xi.
% So we are not going to get the mean medical expenses being exp(m+1/2*sigma^2)
% as footnote 3 of DFJ2010 says, unless stddev(psi)=1. But given the AR(1)
% on zeta and the normal dist of xi,
% stddev(psi)=sqrt(stddev(zeta)^2+stddev(xi)^2)
% [as they are independent normals; AR(1) is asymptotically normal], and then 
% stddev(psi)=sqrt((stddev(zeta innovations)/(1-rho_m))^2+stddev(xi)^2)
% which is clearly not 1.
% So, I introduce a scaling parameter
Params.normalizepsi=sqrt( Params.sigma_mepsilon^2/(1-Params.rho_m^2) + Params.sigma_mxi^2 );
% [Corrected 2026-09-18: the variance of an AR(1) is sigma_eps^2/(1-rho^2), not (sigma_eps/(1-rho))^2. With DFJ2010's
%  numbers this gives 0.050/(1-0.922^2)+0.665=0.999: their psi is unit-variance BY CONSTRUCTION (the C code sets the AR(1)
%  share fracar1 and the white-noise share fracwn=1-fracar1 to sum to one), so normalizepsi is 1.0 and the earlier formula,
%  which gave 2.98, was shrinking the medical expense shock to a third of its size.]
% And then when psi is used anywhere, I replace it with psi/normalizepsi
% Thus stddev(psi/normalizpsi)=1, and we get the mean medical expenses
% being exp(m+1/2*sigma^2).
% [This was not mentioned anywhere in DFJ2010 paper, makes sense as soon as you see it
% because they are two separate estimates off 'same' data, so you cannot just 
% multiply their std devs, but took some figuring out.]




%% Health transition probabilities
% copy-paste from DFJ2010 docs:
% "the parameters of the health transition matrix (healthprof.out)"
% "file describing how to map the two-year transition probabilities that we estimate into one-year transition probabilities (two_to_one.pdf)"

% Start with getting an interpreting data from healthprof.out (copy-paste from DFJ2010 docs):
% "The first column in healthprof.out is age1, which is just an index.
% The remaining columns are logit coefficients in the equation:
% Pr(health_t= bad | X_{t – 2 years}) = exp(X_t b)/(1+exp(X_t b)), 
% The coefficients are as follows: 
% ageshiftFINAL = age-specific intercept term
% healshift = age-specific coefficients for being in bad health 2 years ago (good health is the omitted category)
% maleshift = age-specific coefficients for being male (female is the omitted category)
% PIshift = age-specific coefficients on PI, which is measured as the individual s percentile rank (0 to 1) in the PI distribution 
% bPI2 = coefficients on PI percentile squared.  For now there is no variation across ages, but that can be added.
% Suppose you had an 80-year-old male who was in good health and at the 60th percentile of the PI distribution.
% The probability that this man will be in bad health at age 82 is: 
% Xb = ageshiftFINAL82+maleshift82+PI82*.6+bPI82*.36, 
% where the "82” subscript denotes the row with age index number 82."

% Following is copy-pasted from contents of healthprof.out
healthprof.age=70:1:102; % age index number
healthprof.ageshiftFinal=[-.7599814,-.7719013,-.7743046,-.7680932,-.7541689,-.7334337,-.7067895,-.675138,-.6393813, -.6004211, -.5591595, -.5164981, -.473339, -.4305839, -.3891349, -.3498937, -.3137623, -.2816426, -.2544363, -.2330454, -.2183719, -.2113174, -.2127841, -.2236736, -.2448879, -.277329, -.3218986, -.3794986, -.451031, -.5373976, -.6395003, -.7582409, -.8945215]; 
healthprof.healshift=[2.600759,2.55415,2.507542,2.460934,2.414325,2.367717,2.321109,2.274501,2.227893, 2.181284,2.134676, 2.088068, 2.04146, 1.994851, 1.948243, 1.901635, 1.855027, 1.808418, 1.76181,1.715202,1.668594, 1.621985, 1.575377, 1.528769, 1.482161, 1.435552, 1.388944, 1.342336, 1.295728, 1.249119, 1.202511, 1.155903, 1.109295];
healthprof.maleshift=[.3348944,.3210961,.3072977,.2934994,.2797011,.2659027, .2521044, .2383061, .2245077, .2107094, .1969111, .1831127, .1693144, .1555161, .1417177, .1279194, .1141211, .1003227, .0865244, .0727261, .0589277, .0451294, .0313311, .0175327, .0037344, -.0100639, -.0238623, -.0376606, -.0514589, -.0652573, -.0790556,-.0928539, -.1066523];
healthprof.PIshift=[-1.270953, -1.259267, -1.247581, -1.235895, -1.224209, -1.212523, -1.200836, -1.18915, -1.177464, -1.165778, -1.154092, -1.142406, -1.13072, -1.119033, -1.107347, -1.095661, -1.083975, -1.072289, -1.060603, -1.048917, -1.037231, -1.025545, -1.013858, -1.002172, -.9904861, -.9788, -.9671139, -.9554277, -.9437416, -.9320555, -.9203693, -.9086832, -.8969971];
healthprof.bPI2=.3519772*ones(1,33);
% Note: healthprof.age not use for anything [I just keep it here as it is part of healthprof.out]

% So need to use these to evaluate Xb, then use this to evaluate
% Pr(health=bad|X, two years ago); note gives Prob for 'current age'
probbadhealth_2yr=zeros(2,102-69,2,5); % current h, age (at end of two yr period), gender, income q
for h_c=1:2
    for age=70:102
        for g_c=1:2 % female, male
            for q_c=1:5 % income quintile
                ageindex=age-69;
                cpercentile=(20*(q_c-1)+10)/100; % should be 0.6 is 60th percentile; I identify each quintile with its midpoint, so bottom quintile is 10th percentile [not clear what DFJ2010 did]
                Xb=healthprof.ageshiftFinal(ageindex)+healthprof.healshift(ageindex)*(h_c-1)+healthprof.maleshift(ageindex)*(g_c-1)+healthprof.PIshift(ageindex)*cpercentile+healthprof.bPI2(ageindex)*(cpercentile)^2;
                probbadhealth_2yr(h_c,ageindex,g_c,q_c)=exp(Xb)/(1+exp(Xb)); % Note: exp(x)/(1+exp(x)) is always between 0 and 1 (0 as x--->-Inf, 1 as x--->Inf)
            end
        end
    end
end

% Convert two-year transitions into one-year transitions.
% Idea: interpret two-year transtion from good-to-good health (gg) and two-year 
% transition from bad-to-bad health (bb) as coming from annual (g) and (b). 
% Then a bit of algebra, spelled out in two-to-one.pdf (of
% DJF2010 replication materials) gives us this following.
% Let 
%   A=2-gg-bb
%   B=2*(gg-1)
%   C=1-bb-gg+bb^2
% and then we get that 
%   b=(-B+sqrt(B^2-4AC))/2A  [note, quadratic formula]
% and we can calculate g=(1-bb-b*(1-b))/(1-b)
% As double checks that solution is correct, we should have:
%   gg=g^2+(1-g)*(1-b)

% probbadhealth_2yr=zeros(2,102-69,2,5); % current h, age, gender, income q
probhealth=zeros(2,2,102-69,2,5); % previous h, current h, age, gender, income q
for age=70:102
    for g_c=1:2 % female, male
        for q_c=1:5 % income quintile
            ageindex=age-69;

            hlag_c=1; % h=0 good health
            gg=1-probbadhealth_2yr(hlag_c,ageindex,g_c,q_c);
            hlag_c=2; % h=1 bad health
            bb=probbadhealth_2yr(hlag_c,ageindex,g_c,q_c);

            A=2-gg-bb;
            B=2*(gg-1);
            C=1-bb-gg+bb^2;
            b=(-B+sqrt(B^2-4*A*C))/(2*A);
            % balt=(-B-sqrt(B^2-4AC))/2A % there is another root, but
            % the pdf 'two_to_one.pdf' suggests this is not used
            g=(1-bb-b*(1-b))/(1-b);

            % Do transitions from good health to good/bad health
            hlag_c=1; % h=0 good health
            h_c=1; % h=0 good health
            probhealth(hlag_c,h_c,ageindex,g_c,q_c)=g; % good-to-good
            h_c=2; % h=1 bad health
            probhealth(hlag_c,h_c,ageindex,g_c,q_c)=1-g; % good-to-bad

            % Do transitions from bad health to good/bad health
            hlag_c=2; % h=1 bad health
            h_c=1; % h=0 good health
            probhealth(hlag_c,h_c,ageindex,g_c,q_c)=1-b; % bad-to-good
            h_c=2; % h=1 bad health
            probhealth(hlag_c,h_c,ageindex,g_c,q_c)=b; % bad-to-bad

        end
    end
end
% Note: not entirely clear from DFJ2010 docs which 2yr period becomes which
% 1yr transition (conversion is clear, just not which period we allocate
% the result to). I assume that it is based on the age that we are looking
% at transitions to (not from).

% Done

%% Conditional survival probabilities
% parameters of the survival probability transitions are in deathprof.out

% Following is copy-paste from DFJ2010 docs:
% 2. The line at the bottom of death.do will help you understand deathprof.out
% it works just like healthprof.out: 
% order age1 ageshiftFINAL healshift maleshift PIshift bPI2; 
% The first column in deathprof.out is age1, which is an index.  
% The remaining columns are logit coefficients in the equation:
% Pr(alive at time t |Xt – 2 years, alive at time t – 2 years) = exp(Xtb)/(1+exp(Xtb)), 
% These coefficients are completely analogous to the coefficients of the health 
% transition equation.
% 3. Note that these are 2-year survival probabilities – take the square root of the 2-year survival probability to get the 1-year survival probability.

% Following is copy-pasted from contents of deathprof.out
deathprof.age=70:1:102;
deathprof.ageshiftFinal=[2.546009, 2.565764, 2.575288, 2.574989, 2.565275, 2.546554, 2.519237, 2.483731, 2.440444, 2.389787, 2.332167, 2.267992, 2.197673, 2.121617, 2.040233, 1.95393, 1.863116, 1.7682, 1.669591, 1.567698, 1.462929, 1.355692, 1.246397, 1.135452, 1.023266, .9102477, .7968054, .6833481, .5702842, .4580225, .3469718, .2375407, .1301378];
deathprof.healshift=[-1.026724, -1.012185, -.9976472, -.9831089, -.9685707, -.9540324, -.9394941, -.9249559, -.9104176, -.8958793, -.8813411, -.8668028, -.8522645, -.8377263, -.823188, -.8086497, -.7941114, -.7795732, -.7650349, -.7504966, -.7359584, -.7214201, -.7068818, -.6923436, -.6778053, -.663267, -.6487288, -.6341905, -.6196522, -.605114, -.5905757, -.5760374, -.5614991];
deathprof.maleshift=[-1.017154, -.9981147, -.9790752, -.9600356, -.9409961, -.9219566, -.902917, -.8838775, -.864838, -.8457984, -.8267589, -.8077194, -.7886798, -.7696403, -.7506008, -.7315612, -.7125217, -.6934822, -.6744426, -.6554031, -.6363636, -.6173241, -.5982845, -.579245, -.5602055, -.5411659, -.5221264, -.5030869, -.4840473, -.4650078, -.4459683, -.4269287, -.4078892];
deathprof.PIshift=[.879185, .838859, .7985329, .7582068, .7178808, .6775548, .6372288, .5969027, .5565767, .5162507, .4759246, .4355986, .3952726, .3549465, .3146205, .2742945, .2339684, .1936424, .1533164, .1129903, .0726643, .0323383, -.0079878, -.0483138, -.0886398, -.1289659, -.1692919, -.2096179, -.249944, -.29027, -.3305961, -.3709221, -.4112481];
deathprof.bPI2=.1008181*ones(1,33);

% So need to use these to evaluate Xb, then use this to evaluate
% Pr(alive at t|X, two years ago); note gives Prob for 'current age'
probalive_2yr=zeros(2,102-70+1,2,5); % current h, age, gender, income q
for h_c=1:2
    for age=70:102
        for g_c=1:2 % female, male
            for q_c=1:5 % income quintile
                ageindex=age-69;
                cpercentile=(20*(q_c-1)+10)/100; % should be 0.6 is 60th percentile; I identify each quintile with its midpoint, so bottom quintile is 10th percentile [not clear what DFJ2010 did]
                Xb=deathprof.ageshiftFinal(ageindex)+deathprof.healshift(ageindex)*(h_c-1)+deathprof.maleshift(ageindex)*(g_c-1)+deathprof.PIshift(ageindex)*cpercentile+deathprof.bPI2(ageindex)*(cpercentile);
                probalive_2yr(h_c,ageindex,g_c,q_c)=exp(Xb)/(1+exp(Xb));
            end
        end
    end
end

% Convert to annual by taking the square root
probalive_1yr=sqrt(probalive_2yr); % current h, age, gender, income q


%% Combine health probabilities and conditional survival probabilities to
% construct the transtion matrix for pi_h_J for each ptype.


pi_h_J.femaleq1=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=1; % 1st income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq1(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq1(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq1(3:4,4,:)=1;

% Rest is just repeating this for the other ptypes
pi_h_J.femaleq2=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=2; % 2nd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq2(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq2(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq2(3:4,4,:)=1;

pi_h_J.femaleq3=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=3; % 3rd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq3(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq3(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq3(3:4,4,:)=1;

pi_h_J.femaleq4=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=4; % 4th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq4(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq4(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq4(3:4,4,:)=1;

pi_h_J.femaleq5=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=5; % 5th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq5(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq5(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq5(3:4,4,:)=1;


% Now do male
pi_h_J.maleq1=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=1; % 1st income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq1(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq1(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq1(3:4,4,:)=1;

pi_h_J.maleq2=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=2; % 2nd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq2(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq2(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq2(3:4,4,:)=1;

pi_h_J.maleq3=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=3; % 3rd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq3(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq3(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq3(3:4,4,:)=1;

pi_h_J.maleq4=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=4; % 4th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq4(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq4(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq4(3:4,4,:)=1;

pi_h_J.maleq5=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=5; % 5th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq5(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq5(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq5(3:4,4,:)=1;


%% Income
% (i) the parameters of the income regression (incprof.out)
% The left-hand variable in the income regression is the log of income.  
% The right-hand-side variables are just like those for survival and health:
% age ageshiftFINAL healshift maleshift PIshift bPI2

% Following is copy-pasted from contents of incprof.out
incprof.age=70:1:102;
incprof.ageshiftFinal=[8.024719, 8.031113, 8.035824, 8.039002, 8.040796, 8.041354, 8.040826, 8.039358, 8.037103, 8.034205, 8.030817, 8.027085, 8.023159, 8.019187, 8.015318, 8.011701, 8.008484, 8.005816, 8.003847, 8.002726, 8.002598, 8.003615, 8.005926, 8.009678, 8.01502, 8.022102, 8.031073, 8.042079, 8.055271, 8.070798, 8.088807, 8.109447, 8.132869];
incprof.healshift=0*ones(1,33);
incprof.maleshift=.019004*ones(1,33);
incprof.PIshift=2.458653*ones(1,33);
incprof.bPI2=-.4810466*ones(1,33);

% Create Params from these for use in model
% healshift is 0, so omit from formula in the following (note, female=0, male=1 in the regression of DFJ2010)
cpercentile=10/100; % q1
Params.earnings.femaleq1=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq1=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=30/100; % q2
Params.earnings.femaleq2=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq2=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=50/100; % q3
Params.earnings.femaleq3=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq3=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=70/100; % q4
Params.earnings.femaleq4=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq4=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=90/100; % q5
Params.earnings.femaleq5=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq5=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);


%% Medical expenses: m_coeff and sigma_coeff
% (i) the parameters of the medical expense regression where the left hand 
% side variable is average medical expenses over two years (medexprofX.out); 
% (ii) the parameters of the medical expense regression, where the left hand 
% side variable is annual medical expenses (medexprof_adj.out)

% Copy-paste from DFJ2010 docs
% The list of variables in medexprofX.out and medexprof_adj.out is:
% age ageshiftFINAL healshift maleshift PIshift bPI2 vage vheal vmale vPI vPI2
% This list contains the results from two regressions.  
% The first regression is structured identically to the income regression, 
% except that the left-hand-side variable is log medical expenses.  
% The results of this regression are in the variables
% ageshiftFINAL healshift maleshift PIshift bPI2 
% which are just as before.  
% The second regression is a regression of the variance of medical expenses 
% on the same covariates.  The results of this regression are in the variables
% vage vheal vmale vPI vPI2
% where: 'vage' is an age-specific intercept for the (conditional on X) variance 
% of log medical expenditures; "vheal" is the age-specific set of coefficients 
% for being in bad health; and so on.

% Comment: so first regression is m_coeff, and second regression is the VARIANCE, sigma_coeff=sqrt(second regression)
% Note: the m_coeff is not the mean medical expenses, because of the
% log-expenses formulation the mean expenses is E(m) = exp(mu + 1/2 *sigma^2).

% Following is copy-pasted from contents of medexprof_adj.out
% m_coeff (what DFJ2010 docs refer to as first regression)
medexprof_adj.m_coeff.age=70:1:102;
medexprof_adj.m_coeff.ageshiftFinal=[5.54198301, 5.582088892, 5.609991653, 5.627357155, 5.635732872, 5.636611686, 5.63136711, 5.62131496, 5.607654811, 5.591529413, 5.573962382, 5.555921485, 5.538257273, 5.521761574, 5.507109874, 5.494920062, 5.485692578, 5.479871309, 5.477782758, 5.479694619, 5.485758686, 5.49606959, 5.510602932, 5.529278296, 5.551898496, 5.578208324, 5.607833954, 5.640346602, 5.675199651, 5.711786678, 5.74938749, 5.787221871, 5.824393202];
medexprof_adj.m_coeff.healshift=[-0.126597891, -0.117190485, -0.107783, -0.098375494, -0.088967987, -0.079560502, -0.070152996, -0.06074549, -0.051337984, -0.041930499, -0.032522993, -0.023115486, -0.01370798, -0.004300474, 0.005107011, 0.014514517, 0.023922023, 0.033329508, 0.042737015, 0.052144521, 0.061552006, 0.070959512, 0.080367018, 0.089774503, 0.099182009, 0.108589516, 0.117997001, 0.127404428, 0.136812013, 0.146219435, 0.155626947, 0.165034347, 0.174442059];
medexprof_adj.m_coeff.maleshift=[-0.011614417, -0.021258834, -0.03090325, -0.040547567, -0.050191984, -0.0598363, -0.069480738, -0.079125155, -0.088769371, -0.098413888, -0.108058205, -0.117702642, -0.127346959, -0.136991376, -0.146635692, -0.156280009, -0.165924425, -0.175568763, -0.18521328, -0.194857596, -0.204502013, -0.21414633, -0.223790646, -0.233435084, -0.243079401, -0.252723817, -0.262368134, -0.27201265, -0.281656967, -0.291301405, -0.300945721, -0.310590038, -0.320234455];
medexprof_adj.m_coeff.PIshift=[1.707721202, 1.794494695, 1.881268188, 1.968040892, 2.054814385, 2.141587878, 2.228361583, 2.315135076, 2.401907781, 2.488681274, 2.575454766, 2.662227471, 2.749000964, 2.835774457, 2.922547162, 3.009320655, 3.096093147, 3.182866852, 3.269640345, 3.35641405, 3.443187543, 3.529960036, 3.61673374, 3.703507233, 3.790280726, 3.877053431, 3.963826924, 4.050599629, 4.137373122, 4.224146614, 4.310920319, 4.397692812, 4.484466305];
medexprof_adj.m_coeff.bPI2=-1.537435863*ones(1,33);
% sigma_coeff (what DFJ2010 docs refer to as second regression)
medexprof_adj.sigma_coeff.age=70:1:102;
medexprof_adj.sigma_coeff.ageshiftFinal=[0.859672381, 0.904219617, 0.934751994, 0.95437459, 0.965906155, 0.972021629, 0.97510948, 0.977415081, 0.980896477, 0.987369474, 0.998361836, 1.01525923, 1.039159854, 1.071018951, 1.111505152, 1.161143277, 1.220172244, 1.288686881, 1.366495683, 1.453264763, 1.548370627, 1.651049821, 1.760250136, 1.874778407, 1.993151008, 2.113743353, 2.234646092, 2.353801797, 2.468869698, 2.577366643, 2.67651902, 2.763412258, 2.834845595];
medexprof_adj.sigma_coeff.healshift=[0.641515082, 0.66864537, 0.6957758, 0.722906087, 0.750036375, 0.777166805, 0.804297092, 0.83142738, 0.858557668, 0.885688098, 0.912818385, 0.939948673, 0.96707896, 0.994209248, 1.021339678, 1.048469965, 1.075600253, 1.102730683, 1.129860971, 1.156991258, 1.184121688, 1.211251976, 1.238382263, 1.265512693, 1.292642981, 1.319773269, 1.346903699, 1.374033844, 1.401164274, 1.428295131, 1.455424707, 1.482555706, 1.509685282];
medexprof_adj.sigma_coeff.maleshift=[0.329451635, 0.307158368, 0.284865101, 0.262571834, 0.240278567, 0.217985301, 0.195692176, 0.173398909, 0.151105643, 0.128812376, 0.106519109, 0.084225985, 0.061932718, 0.039639451, 0.017346184, -0.004947083, -0.027240349, -0.049533474, -0.071826741, -0.094120007, -0.116413274, -0.138706541, -0.160999808, -0.183292932, -0.205586199, -0.227879466, -0.250172732, -0.272465999, -0.294759266, -0.31705239, -0.339345657, -0.361638924, -0.383932191];
medexprof_adj.sigma_coeff.PIshift=[1.460772596, 1.51177261, 1.562772625, 1.613771215, 1.664771229, 1.715771244, 1.766769834, 1.817769848, 1.868768439, 1.919768453, 1.970768467, 2.021767058, 2.072767072, 2.123767086, 2.174765677, 2.225765691, 2.276765705, 2.327764295, 2.37876431, 2.4297629, 2.480762914, 2.531762929, 2.582761519, 2.633761533, 2.684761548, 2.735760138, 2.786760152, 2.837758743, 2.888758757, 2.939758771, 2.990757362, 3.041757376, 3.09275739];
medexprof_adj.sigma_coeff.bPI2=-1.384958574*ones(1,33);

% Create Params from these for use in model
% Note, we need two separate Params for health, and can use if-statement in
% ReturnFn to figure out which one applies for current health state
% note: female=0 and male=1 for maleshift (following DFJ2010 regressions)
% note: good health=0 and bad health =1 for healshift (following DFJ2010 regressions)
cpercentile=10/100; % q1
Params.m_coeff_healthgood.femaleq1=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq1 =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq1  =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq1   =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=30/100; % q2
Params.m_coeff_healthgood.femaleq2=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq2 =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq2  =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq2   =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=50/100; % q3
Params.m_coeff_healthgood.femaleq3=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq3 =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq3  =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq3   =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=70/100; % q4
Params.m_coeff_healthgood.femaleq4=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq4 =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq4  =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq4   =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=90/100; % q5
Params.m_coeff_healthgood.femaleq5=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq5 =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq5  =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq5   =medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;



cpercentile=10/100; % q1
Params.sigma_coeff_healthgood.femaleq1=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq1 =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq1  =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq1   =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=30/100; % q2
Params.sigma_coeff_healthgood.femaleq2=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq2 =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq2  =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq2   =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=50/100; % q3
Params.sigma_coeff_healthgood.femaleq3=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq3 =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq3  =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq3   =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=70/100; % q4
Params.sigma_coeff_healthgood.femaleq4=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq4 =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq4  =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq4   =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=90/100; % q5
Params.sigma_coeff_healthgood.femaleq5=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq5 =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq5  =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq5   =medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
% The second regression predicts the VARIANCE of log medical expenses (DFJ2010 docs: 'vage is an age-specific intercept for
% the (conditional on X) variance of log medical expenditures'), and the GAUSS code takes its square root before use
% (getdat25_wave8.g getmxtab: stdmxs=sqrt(mxvarscl*stdmxs); wealth20_Singles.c line 1627: sqrt(hcSig^2+PIC41*PI1+PIC42*PI2)).
% medexprof_adj.out already has the mxvarscl (two-year to one-year) factor applied, so the standard deviation is the square
% root of the prediction. [Corrected 2026-09-18: earlier versions used the variance itself as the standard deviation, which
% overstates medical expense risk at old ages (a variance of 3 at age 95 was being used as a std dev of 3 rather than 1.7).]
for ii=1:length(Names_i)
    Params.sigma_coeff_healthgood.(Names_i{ii})=sqrt(Params.sigma_coeff_healthgood.(Names_i{ii}));
    Params.sigma_coeff_healthbad.(Names_i{ii})=sqrt(Params.sigma_coeff_healthbad.(Names_i{ii}));
end


