function dm = DIAMOBIL_TEST(dg, da)
% "DIAMOBIL" calculates the mobility diameter by interpolation between...
%   ...free molecular and continuum limit approximations.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     dg: Aggregates' gyration diameter
%     da: ~ projected area equivalent diameter
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     dm: Mobility diameter
% ----------------------------------------------------------------------- %

% constant coefficients for slip correction
A1 = 1.165;
A2 = 0.483;
A3 = 0.997;

lambda = 6.8e-8; % mean free path

mu = 18.37e-6; % dynamic viscosity

% Cunningham correction factor (x: particle diameter)
cc = @(x) 1 + x .* (A1 + A2 * exp(-A3 ./ x));

dm_c = 0.75 * dg; % continuum regime approximation initialized

% dm_c = zeros(length(dg),1);
% dm_c(npp < 100) = dpp(npp < 100) .* npp(npp < 100).^0.46;
% dm_c(npp >= 100) = 0.65 * dpp(npp >= 100) .* npp(npp >= 100).^0.56;

R_adj = da.^2 ./ (2 * dm_c);

f_not = (mu * 3 * pi * dm_c) ./ cc(lambda ./ R_adj);

dm0 = 0.75 * dg; % initial guess

dm = zeros(length(dm0),1);
dm2 = zeros(length(dm0),1);

figure(1)

for i = 1 : length(dm0)
    f = @(x) 3 * pi * mu * x ./ cc(2 * lambda ./ x);
    df = @(x) 3 * pi * mu * x ./ cc(2 * lambda ./ x) - f_not(i);
    
    dm(i) = fzero(df,(dm0(i) + da(i)) / 2);
%     dm(i) = lsqnonlin(df,(dm0(i) + da(i)) / 2);
    
    dm_fit = [0.1 * min(dm_c(i),da(i)), 10 * max(dm_c(i),da(i))];
    dm_fit = logspace(log10(dm_fit(1)), log10(dm_fit(2)), 1e4);
    
    f_fit = f(dm_fit);
    df_fit = df(dm_fit);
    dm2(i) = dm_fit(abs(df_fit) == min(abs(df_fit)));
    
    plot(dm_fit, f_fit)
    hold on
    plot(dm_fit, repmat(f_not(i), 1e4, 1))
    hold off
    legend('f_{d_m}','f_{d_{adj}}')
    
    set(gca, 'XScale', 'log', 'YScale', 'log')
    axis([-inf inf -inf inf])
    xlabel('d_m')
    ylabel('f')
end

end

