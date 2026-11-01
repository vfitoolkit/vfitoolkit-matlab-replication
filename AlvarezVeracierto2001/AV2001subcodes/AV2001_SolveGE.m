function [i_eqm,Out]=AV2001_SolveGE(ilo,ihi,Params,Grids,vfoptions,simoptions,verbose,fzerooptions_i,label)
% Solves for the equilibrium interest rate: brackets the goods-market residual
% of AV2001_StationaryEqm on [ilo,ihi] (widening if needed), then fzero.
% label is used in messages.

reslo=AV2001_StationaryEqm(ilo,Params,Grids,vfoptions,simoptions,verbose);
reshi=AV2001_StationaryEqm(ihi,Params,Grids,vfoptions,simoptions,verbose);
bracketiter=0;
while reslo*reshi>0 && bracketiter<10
    ilo=max(ilo-0.001,0.0005); ihi=ihi+0.001;
    reslo=AV2001_StationaryEqm(ilo,Params,Grids,vfoptions,simoptions,verbose);
    reshi=AV2001_StationaryEqm(ihi,Params,Grids,vfoptions,simoptions,verbose);
    bracketiter=bracketiter+1;
end
if reslo*reshi>0
    error('AV2001_SolveGE (%s): could not bracket the equilibrium interest rate',label)
end
i_eqm=fzero(@(ii) AV2001_StationaryEqm(ii,Params,Grids,vfoptions,simoptions,verbose),[ilo,ihi],fzerooptions_i);
[~,Out]=AV2001_StationaryEqm(i_eqm,Params,Grids,vfoptions,simoptions,verbose);

end
