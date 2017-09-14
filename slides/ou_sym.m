%
% GPSS'17 appendix -
% Symbolic computations for the covariance function and
% transition density of the OU process. By SS'17. 
%

    %%
    % Compute the covariance function
    %
    syms lam q positive;
    syms tau w real;
    
    Sw = q/(lam^2 + w^2);
    
    disp('S(w) = ');
    pretty(Sw)
    
    disp('K(tau) = ');
    Ktau = ifourier(Sw,tau);
    pretty(Ktau)
    
    %%
    % Compute the transition density. The density is then
    % p(f(t + Dt) | f(t)) = N(f(t + Dt) | A f(t), Q).
    %
    syms dt q lam s positive;
    
    A = exp(-lam*dt)
    Q = int(exp(-2*lam*(dt-s))*q,0,dt)
    
    
    
    