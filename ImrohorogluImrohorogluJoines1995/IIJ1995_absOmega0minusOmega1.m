function EqmCondns=IIJ1995_absOmega0minusOmega1(GEandLumpSum,Omega1,Params,jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, FnsToEvaluate2, GeneralEqmEqns, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, FnsToEvaluateParamNames2, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% Welfare with compensation L in each period of life and without social security

Params.r=GEandLumpSum(1);
Params.Tr_beq=GEandLumpSum(2);
Params.LumpSum=GEandLumpSum(3);

GEprices=GEandLumpSum(1:2);

N_z=prod(n_z);

if heteroagentoptions.verbose==1
    fprintf('Currently evaluating: LumpSum=%8.4f \n', Params.LumpSum)
end


%% 
% If 'exogenous shock fn' is used and depends on GE parameters then
% precompute it here (otherwise it is already precomputed).
if isfield(vfoptions,'ExogShockFn')
    if ~isfield(vfoptions,'pi_z_J') % This is implicitly checking that ExogShockFn does depend on GE params (if it doesn't then this field will already exist)
        pi_z_J=zeros(N_z,N_z,N_j);
        for jj=1:N_j
            if isfield(vfoptions,'ExogShockFnParamNames')
                ExogShockFnParamsVec=CreateVectorFromParams(Parameters, simoptions.ExogShockFnParamNames,jj);
                ExogShockFnParamsCell=cell(length(ExogShockFnParamsVec),1);
                for ii=1:length(ExogShockFnParamsVec)
                    ExogShockFnParamsCell(ii,1)={ExogShockFnParamsVec(ii)};
                end
                [z_grid,pi_z]=simoptions.ExogShockFn(ExogShockFnParamsCell{:});
            else
                [z_grid,pi_z]=simoptions.ExogShockFn(jj);
            end
            pi_z_J(:,:,jj)=gather(pi_z);
            z_grid_J(:,jj)=gather(z_grid);
        end
        % Now store them in vfoptions and simoptions
        vfoptions.pi_z_J=pi_z_J;
        vfoptions.z_grid_J=z_grid_J;
        simoptions.pi_z_J=pi_z_J;
        simoptions.z_grid_J=z_grid_J;
    end
end

%% 


[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

% Calculate the general eqm conditions
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);
GeneralEqmConditionsVec=real(GeneralEqmConditions_Case1(AggVars,GEprices, GeneralEqmEqns, Params,GeneralEqmEqnParamNames));

% Calculate absOmega0minusOmega1
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate2, Params, FnsToEvaluateParamNames2, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid,[],simoptions);
UtilityOnGrid=shiftdim(ValuesOnGrid(5,:,:,:),1);
discountongrid=shiftdim(cumprod(Params.beta*Params.sj),-1);
AgeConditionalStationaryDist=StationaryDist./sum(sum(StationaryDist,1),2);
Omega0=sum(sum(sum(discountongrid.*AgeConditionalStationaryDist.*UtilityOnGrid)));
absOmega0minusOmega1=abs(Omega0-Omega1);

EqmCondns=[gather(abs(GeneralEqmConditionsVec)),gather(absOmega0minusOmega1)];

if heteroagentoptions.verbose==1
    fprintf('GEprices: r=%8.4f Tr_beq=%8.4f \n', Params.r, Params.Tr_beq)
    fprintf('GeneralEqmConditions: %8.4f %8.4f \n', GeneralEqmConditionsVec)
    fprintf('Omega0 and Omega1: Omega0=%8.4f Omega1=%8.4f \n', Omega0, Omega1)
end

EqmCondns=sum([1,10,1].*(EqmCondns)); % I put a larger weight on Tr_beq as this is almost a perfect substitute for LumpSum



end
