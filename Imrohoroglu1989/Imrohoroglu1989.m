% Imrohoroglu (1989) - Cost of Business Cycles with Indivisibilities and Liquidity Constraints

% INCOMPLETE, DOES NOT YET REPLICATE TABLE 1

% You can change this grid size to see if anything changes
npoints_basic_a_grid=301; % Page 1373 of Imrohoroglu (1989), uses 301 (codes assume it is an odd number)


%%
% Create some matrices to store numbers needed for tables
AvgUtility=nan(2,6); % WhichSigma-by-EconomyEnvironment
AvgConsumption=nan(2,6);


%Choose one of the following
%EconomyEnvironment pairs
%1: Economy=1, Environment=A
%2: Economy=2, Environment=A
%3: Economy=1, Environment=B
%4: Economy=2, Environment=B
%5: Economy=1, Environment=C
%6: Economy=2, Environment=C
%EconomyEnvironment=1
% Which of the two values of sigma, the risk aversion parameter, to use.
%WhichSigma=1

for WhichSigma=1:2
    for EconomyEnvironment=1:6
        
        
        %% Setup
        
        vfoptions.parallel=2; % Parallelize on GPU (this is anyway the default)
        
        % num_z_vars=1; %This is 'n' in Imrohoroglu
        % num_s_vars=1; %This is 'i' in Imrohoroglu
        % num_a_vars=1; %This is 'a' in Imrohoroglu
        % num_d_vars=0; %There is none in Imrohorglu (you only choose aprime, which Case 1 codes assume is chosen anyway)
        
        %Discounting rate
        Params.beta = 0.995
        
        %Parameters
        if WhichSigma==1
            Params.sigma=1.5;
        elseif WhichSigma==2
            Params.sigma=6.2; %Used for comparison of results with Lucas (1987) which uses 6.2
        end
        Params.theta=0.25;
        Params.r_l=0; %Interest rate on lending (in environment A, this is the only r)  (pg. 1370)
        Params.r_b=0.08; %Interest rate on borrowing (pg. 1370)
        Params.y=1; %Normalization (THIS IS NOT FROM PAPER, I JUST GUESSED)
        
        s_grid=[1,2]';
        if EconomyEnvironment==1 || EconomyEnvironment==2 || EconomyEnvironment==3 % Economy with aggregate shocks
            z_grid=[1,2]';
            sz_grid=[s_grid;z_grid];
            pi_sz=[0.9141,0.0234,0.0587, 0.0038; 0.5625, 0.3750, 0.0269, 0.0356; 0.0608, 0.0016, 0.8813, 0.0563; 0.0375, 0.0250, 0.4031, 0.5344]; % From Eqn 14 of Imrohoroglu (1989)
        elseif EconomyEnvironment==4 || EconomyEnvironment==5 || EconomyEnvironment==6 % Economy without aggregate shocks
            z_grid=[1]';
            sz_grid=[s_grid;z_grid];
            pi_sz=[0.9565,0.0435;0.5,0.5]; % From Eqn 15 of Imrohoroglu (1989)
        end
        
        if EconomyEnvironment==1 || EconomyEnvironment==4 % Environment A
            a_grid=linspace(0,8,npoints_basic_a_grid)'; % Page 1373 of Imrohoroglu (1989)
        elseif EconomyEnvironment==2 || EconomyEnvironment==5 % Environment B
            a_grid=linspace(-8,8,2*npoints_basic_a_grid-1)';
        elseif EconomyEnvironment==3 || EconomyEnvironment==6 % Environment C
            a_grid=0;
        end
        
        n_s=length(s_grid);
        n_z=length(z_grid);
        n_sz=[n_s,n_z];
        n_a=length(a_grid);
        
        %% Now, create the return function
        DiscountFactorParamNames={'beta'};
        
        ReturnFn=@(aprime_val, a_val, s_val,z_val,r_l,r_b,y,theta,sigma) Imrohoroglu1989_ReturnFn(aprime_val, a_val, s_val,z_val,r_l,r_b,y,theta,sigma);
        ReturnFnParamNames={'r_l','r_b','y','theta','sigma'}; %It is important that these are in same order as they appear in 'Imrohoroglu1989_ReturnFn'
        
        
        %% Solve
        %Do the value function iteration. Returns both the value function itself,
        %and the optimal policy function.
        d_grid=0; %no d variable
        n_d=0; %no d variable
        
        tic;
        V0=ones([n_a,n_sz]);
        [V, Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_sz,d_grid,a_grid,sz_grid, pi_sz, DiscountFactorParamNames, ReturnFn, vfoptions, Params, ReturnFnParamNames);
        time=toc;
        
        fprintf('Time to solve the value function iteration was %8.2f seconds. \n', time)
        
        
        %% Draw a graph of the value function
        
        % To get 'nicer' x and y axes use
        if n_z==1 && n_a>1
            surf(a_grid*ones(1,prod(n_s)),ones(n_a,1)*s_grid',V(:,:))            
        elseif n_z==2 && n_a>1
            surf(a_grid*ones(1,prod(n_s)),ones(n_a,1)*s_grid',V(:,:,1))
            hold on
            surf(a_grid*ones(1,prod(n_s)),ones(n_a,1)*s_grid',V(:,:,2))
            hold off
            % It is not so visually obvious but there are two surfaces being plotted.
            % One for z=1 and one for z=2 (they sit almost immediately one on top of
            % the other).
        end
        
        %%
        % Imrohoroglu (1989) does 500000 period simulations (apparently without
        % burnin). Here we instead use the more robust method of iterating on whole distribution.
        StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_sz,pi_sz);
        
        %% Generate some output following what is reported in Imrohoroglu (1989)
        nsample=10^6; npts=301;
        % Figure 1: Asset holdings (Imrohoroglu gives this figure only for
        % Environment A) (For EconomyEnvironment's 3 & 6 there is nothing
        % to plot)
        fig1=figure(1);
        if EconomyEnvironment==1 || EconomyEnvironment==2
            X21=randsample(a_grid,nsample,true,StationaryDist(:,2,1));
            [freq11,X11i] = ksdensity(a_grid,'weights',gather(StationaryDist(:,1,1)/sum(StationaryDist(:,1,1))),'npoints',npts,'kernel','epanechnikov','bandwidth',0.3); %kernel-smooth estimate of pdf
            [freq21,X21i] = ksdensity(a_grid,'weights',gather(StationaryDist(:,2,1)/sum(StationaryDist(:,2,1))),'npoints',npts,'kernel','epanechnikov','bandwidth',0.3);
            subplot(2,1,1); plot(X11i,freq11*sum(StationaryDist(:,1,1)),X21i,freq21*sum(StationaryDist(:,2,1)))
            title('Good Times','FontSize',18);
            
            X12=randsample(a_grid,nsample,true,StationaryDist(:,1,2));
            X22=randsample(a_grid,nsample,true,StationaryDist(:,2,2));
            [freq12,X12i] = ksdensity(a_grid,'weights',gather(StationaryDist(:,1,2)/sum(StationaryDist(:,1,2))),'npoints',npts,'kernel','epanechnikov','bandwidth',0.3);
            [freq22,X22i] = ksdensity(a_grid,'weights',gather(StationaryDist(:,2,2)/sum(StationaryDist(:,2,2))),'npoints',npts,'kernel','epanechnikov','bandwidth',0.3);
            subplot(2,1,2); plot(X12i,freq12*sum(StationaryDist(:,1,2)),X22i,freq22*sum(StationaryDist(:,2,2)))
            title('Bad Times','FontSize',18);
            
            % % A little comparison of which method (out of three different tries) was doing best in terms of graphing
            % % the pdf. ksdensity was clear winner (based on judging against the cdf)
            % % Unfortunately Matlab ksdensity() does not yet appear to allow for data-based automatic bandwidth selection.
                figure(3);
                X11=randsample(a_grid,nsample,true,StationaryDist(:,1,1)); % create random sample X11 based on StationaryDist
                [N11,edges11] = histcounts(X11,npts);
                plot(X11i,freq11*sum(StationaryDist(:,1,1)),edges11(2:end)+edges11(1:end-1)/2,(N11/sum(N11))*sum(StationaryDist(:,1,1)),a_grid,StationaryDist(:,1,1),a_grid,cumsum(StationaryDist(:,1,1)));
        elseif EconomyEnvironment==4 || EconomyEnvironment==5
            plot(a_grid,StationaryDist(:,1),a_grid,StationaryDist(:,2))
        end
        
        if EconomyEnvironment==1 && WhichSigma==1
            set(fig1, 'Color', 'white');     % white bckgr
            set(fig1, 'Unit', 'centimeters');  % set unit to cm
            set(fig1,'position',[0 0 20 10]);  % set size
            export_fig(fig1, ...            % figure handle
                ['./SavedOutput/Graphs/Imrohoroglu1989_Figure1.pdf'],... % name of output file without extension   %  '-painters', ...            % renderer
                '-pdf', ...                 % file format
                '-r300' );                  % resolution
        end
        
        %
        fig2=figure(2);
        PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_sz,d_grid,a_grid,vfoptions.parallel);
        if EconomyEnvironment==1 || EconomyEnvironment==2 || EconomyEnvironment==3
            plot(a_grid,PolicyValues(1,:,1,1)-a_grid',a_grid,PolicyValues(1,:,1,2)-a_grid',a_grid,PolicyValues(1,:,2,1)-a_grid',a_grid,PolicyValues(1,:,2,2)-a_grid')
        elseif EconomyEnvironment==4 || EconomyEnvironment==5 || EconomyEnvironment==6
            plot(a_grid,PolicyValues(1,:,1)-a_grid',a_grid,PolicyValues(1,:,2)-a_grid')
        end
        
        %%
        SSvalueParamNames={}; % Only needed if using gpu for SSvalues
        SSvalue_Utility = @(aprime_val,a_val,sz_val,pi_sz,p_val) Imrohoroglu1989_ReturnFn(aprime_val, a_val, sz_val(1),sz_val(2),Params.r_l,Params.r_b,Params.y,Params.theta,Params.sigma);
        SSvalue_Consumption = @(aprime_val,a_val,sz_val,pi_sz,p_val) Imrohoroglu1989_ConsFn(aprime_val, a_val, sz_val(1),sz_val(2),Params.r_l,Params.r_b,Params.y,Params.theta,Params.sigma);
        SSvalue_AssetsBorrowed = @(aprime_val,a_val,sz_val,pi_sz,p_val,param) a_val*(a_val<0);
        SSvalue_AssetsStored = @(aprime_val,a_val,sz_val,pi_sz,p_val,param) a_val;
        SSvalue_AssetsSaved = @(aprime_val,a_val,sz_val,pi_sz,p_val,param) a_val*(a_val>0);
        SSvalue_Earnings = @(aprime_val,a_val,sz_val,pi_sz,p_val,param) Params.y*((sz_val(1)==1)+Params.theta*(sz_val(1)==2));
        SSvalue_Income = @(aprime_val,a_val,sz_val,pi_sz,p_val,param) Imrohoroglu1989_IncomeFn(aprime_val, a_val, sz_val(1),sz_val(2),Params.r_l,Params.r_b,Params.y,Params.theta,Params.sigma);
        SSvaluesFn={SSvalue_Utility, SSvalue_Consumption, SSvalue_AssetsBorrowed, SSvalue_AssetsStored, SSvalue_AssetsSaved,SSvalue_Earnings,SSvalue_Income};
        
        SSvalues_AggVars=SSvalues_AggVars_Case1(gather(StationaryDist), gather(Policy), SSvaluesFn,Params, SSvalueParamNames,n_d, n_a, n_sz, d_grid, a_grid,sz_grid,pi_sz,[], 1)
        
        AvgUtility(WhichSigma,EconomyEnvironment)=SSvalues_AggVars(1);
        AvgConsumption(WhichSigma,EconomyEnvironment)=SSvalues_AggVars(2);
               
        if WhichSigma==1
            if EconomyEnvironment==1
                Table2(1:5,2)=[SSvalues_AggVars(3),SSvalues_AggVars(4),SSvalues_AggVars(5),SSvalues_AggVars(7),SSvalues_AggVars(2)];%[SSvalue_AssetsBorrowed; SSvalue_AssetsStored; SSvalue_AssetsSaved; SSvalue_Income; SSvalue_Consumption];
            elseif EconomyEnvironment==2
                Table2(1:5,1)=[SSvalues_AggVars(3),SSvalues_AggVars(4),SSvalues_AggVars(5),SSvalues_AggVars(7),SSvalues_AggVars(2)];
            end
        end
        
        
        % End the for loops for WhichSigma and EconomyEnvironment
    end
end


%% Create Table 1
%Table1=nan(2,2);
% NOT YET IMPLEMENTED








