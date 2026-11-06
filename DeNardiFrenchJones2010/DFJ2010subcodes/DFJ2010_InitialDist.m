% DFJ2010_InitialDist: the 1996 joint distribution of the state variables, on the model's grids
%
% DFJ2010 (paper Sec. IV.A, GAUSS initdist()) start every simulated individual from a draw of
% (age, sex, permanent income, health, assets) from the 1996 AHEAD data, restricted to those with
% assets, medical expenses and health all observed in 1996 and with 1996 age in [70,102]. The
% persistent medical-expense shock zeta is drawn from its stationary distribution and xi is iid,
% so those are not data objects (the model script builds them from pi_zeta and pi_xi).
%
% Script operating on workspace variables (set by DFJ2010_DataWork):
%   person (table with female, PI), PIq_in, A96, H96, AGE96, got96 (logical: usable 1996 obs)
%   a_grid, n_a (must match the model script), Names_i
% Produces:
%   InitDist.(ptype)   n_a x 2 x 33 mass over (asset grid point, health {good,bad}, 1996 age 70..102),
%                      summing to one within each ptype; assets are placed on a_grid by a lottery
%                      between the two bracketing grid points (preserves the mean; negative net
%                      worth goes to a=0)
%   AgeDist.(ptype)    1 x 33 distribution of 1996 age within the ptype
%   ptypemasses        1 x 10 share of the 1996 sample in each ptype (Names_i order)
%   InitCounts         1 x 10 number of 1996 observations behind each ptype
%   InitObs            table of the raw 1996 draws (HHIDPN, ptype, age, badhealth, assets)

N_j=33; % ages 70..102
nP=numel(Names_i);
InitDist=struct();
AgeDist=struct();
InitCounts=zeros(1,nP);
InitObs=table();

for pp=1:nP
    name=Names_i{pp};
    isfemale=strncmp(name,'female',6);
    q=str2double(name(end));
    sub=got96 & (person.female==isfemale) & (PIq_in==q);
    n=sum(sub);
    InitCounts(pp)=n;
    D=zeros(n_a,2,N_j);
    ages=AGE96(sub);
    hs=H96(sub);
    as=A96(sub);
    for ii=1:n
        jj=ages(ii)-69; % age 70 -> j=1
        hh=hs(ii)+1;    % good health (0) -> 1, bad health (1) -> 2
        a=as(ii);
        if a<=a_grid(1)
            D(1,hh,jj)=D(1,hh,jj)+1;
        elseif a>=a_grid(end)
            D(end,hh,jj)=D(end,hh,jj)+1;
        else
            il=find(a_grid<=a,1,'last');
            wl=(a_grid(il+1)-a)/(a_grid(il+1)-a_grid(il));
            D(il,hh,jj)=D(il,hh,jj)+wl;
            D(il+1,hh,jj)=D(il+1,hh,jj)+(1-wl);
        end
    end
    if n>0
        D=D/n;
    end
    InitDist.(name)=D;
    AgeDist.(name)=reshape(sum(D,[1,2]),1,N_j);
    InitObs=[InitObs; table(person.HHIDPN(sub),repmat(string(name),n,1),ages,hs,as,'VariableNames',{'HHIDPN','ptype','age1996','badhealth','assets'})]; %#ok<AGROW>
end
ptypemasses=InitCounts/sum(InitCounts);

fprintf('DFJ2010_InitialDist: %d usable 1996 observations (DFJ initdist requires assets, medex and health observed)\n',sum(InitCounts))
fprintf('   %-10s %5s %8s %10s %9s %10s %8s\n','ptype','n','share','medassets','frac a<0','fracbadhlth','meanage')
for pp=1:nP
    name=Names_i{pp};
    sub=strcmp(InitObs.ptype,name);
    fprintf('   %-10s %5d %8.3f %10.0f %9.3f %10.3f %8.1f\n',name,InitCounts(pp),ptypemasses(pp),median(InitObs.assets(sub)),mean(InitObs.assets(sub)<0),mean(InitObs.badhealth(sub)),mean(InitObs.age1996(sub)))
end
fprintf('   mass on a_grid(1)=0 by ptype: %s\n',mat2str(cellfun(@(nm) sum(InitDist.(nm)(1,:,:),'all'),Names_i),3))
