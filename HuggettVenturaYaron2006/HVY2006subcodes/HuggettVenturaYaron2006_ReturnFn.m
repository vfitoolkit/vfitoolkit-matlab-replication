function F=HuggettVenturaYaron2006_ReturnFn(l,h,w)

F=-Inf;

% f_hdl=ability*(h*l)^alpha; % new production of human capital
% hprime=f_hdl+h*(1-delta);

if l>=0 && l<=1 % l has to be between 0 and 1 [l=0 is no investment, which Ben-Porath agents choose near the end of working life]
    F=w*h*(1-l);
end

end