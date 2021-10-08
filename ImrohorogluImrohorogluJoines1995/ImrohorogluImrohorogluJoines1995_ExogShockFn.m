function [z_grid,pi_z]=ImrohorogluImrohorogluJoines1995_ExogShockFn(agej,Jr, MedicalShock)
% See IIJ1995, pg 111

% if I_j==1
z_grid=[0; 1]; % Indicates 'employed'
pi_z=[0.06, 0.94; 0.06, 0.94]; % Note that IIJ1995 reports this on pg 92 with the 'reversed' definition of which state is which.
if agej==Jr
    % When they move from Jr into the next period, which is retirement, we
    % want to redistribute the z shock to the relevant new distribution.
    if MedicalShock==1
        pi_z=[1,0;1,0]; % Everyone starts healthy
        % Originally I used the unconditional probability of illness, but this delivered results at stark odds with paper.
%         pi_z=[0.82,0.18; 0.82,0.18]; % iid Probability of illness is 0.18
    elseif MedicalShock==2
        pi_z=[1,0;1,0]; % Everyone starts healthy
        % Originally I used the unconditional probability of illness, but this delivered results at stark odds with paper.
%         pi_z=[0.91,0.09; 0.91,0.09]; % iid Probability of illness is 0.09
    end
elseif agej>Jr % z now indicates 'illness', or more precisely 'cost of illness as percent of employed wage.
    if MedicalShock==1
        z_grid=[0;0.25]; % Cost of illness is 25 percent of employed wage
        pi_z=[0.9450, 0.0550; 0.25, 0.75]; % Implies stationarty distribution of [0.8197; 0.1903]
    elseif MedicalShock==2
        z_grid=[0;0.35]; % 35 percent of employed wage
        pi_z=[0.9753, 0.0247; 0.25, 0.75]; % Implies stationarty distribution of [0.9101; 0.0899]
    end
end

end