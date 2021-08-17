function GeneralEqmConditions_and_absOmega0minusOmega1=IIJ1995_absOmega0minusOmega1(GEprices_and_L,Omega1,Params,jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, FnsToEvaluate2, GeneralEqmEqns, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, FnsToEvaluateParamNames2, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% Welfare with compensation L in each period of life and without social security

Params.b=GEprices_and_L(1);
Params.tau_u=GEprices_and_L(2);
Params.tau_s=GEprices_and_L(3);
Params.Tr_beq=GEprices_and_L(4);
Params.LumpSum=GEprices_and_L(5);

GEprices=GEprices_and_L(1:4);

N_z=prod(n_z);

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

% Get things ready to output and print progress
GeneralEqmConditionsVec=[GeneralEqmConditionsVec,absOmega0minusOmega1];

heteroagentoptions.multiGEweights=[1,1,1,1,5]; % The magnitude of absOmega0minusOmega1 is small so I multiply it by 5 to ensure it is given some importance.
GeneralEqmConditions_and_absOmega0minusOmega1=sum(abs(heteroagentoptions.multiGEweights.*GeneralEqmConditionsVec));

GeneralEqmConditions_and_absOmega0minusOmega1=gather(GeneralEqmConditions_and_absOmega0minusOmega1);

if heteroagentoptions.verbose==1
    fprintf('Current GE prices and GeneralEqmConditionsVec, also LumpSum and absOmega0minusOmega1 respectively. \n')
    GEprices_and_L
    GeneralEqmConditionsVec
end





end
