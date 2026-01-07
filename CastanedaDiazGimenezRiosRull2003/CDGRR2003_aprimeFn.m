function aprime=CDGRR2003_aprimeFn(d,z,zprime,J,tauE,zlowerbar)
% interitanceasset: a'(d,z,z')

% If you don't die, then you just choose aprime directly
aprime=d;
% If you die, then you get what is left after the estate tax
if z>J
    if d>zlowerbar
        aprime=zlowerbar+(1-tauE)*(d-zlowerbar);
    % else: just keep aprime=d
    end
end


end