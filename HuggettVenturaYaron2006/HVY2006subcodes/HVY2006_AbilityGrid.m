function Params=HVY2006_AbilityGrid(Params)
% ability and abilitydist are parametrized by mean_logability and stddev_logability

[ability_grid,pi_ability]=discretizeAR1_FarmerToda(Params.mean_logability,0,Params.stddev_logability,Params.Number_i);
ability_grid=exp(ability_grid);
pi_ability=pi_ability(1,:)'; % iid
% Comment: originally I just created a grid on ability, rather than log of
% ability, but this gave some points a negative value for ability, so
% switched to current setup of creating the grid on log of ability, then
% take exponent.
% Note that the initial distribution is anyway log-normal, so working with
% a grid on log-ability makes more sense for that anyway.

% Ability is a parameter that depends on permanent type
Params.ability=ability_grid;

% PTypeDistParamNames={'abilitydist'};
Params.abilitydist=pi_ability;


end