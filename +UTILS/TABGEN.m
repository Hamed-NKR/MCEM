function tab = TABGEN(dpp, dm)
% "TABGEN" generates a table of properties for typical fractal soot...
%   ...according to the well-established literature correlations that...
%   ...can be found in Sipkens et al. (2023)'s review on soot morphology.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   dpp: Mean primary particle diameter
%   dm: Mobility diameter
% ----------------------------------------------------------------------- %
%
% Outputs:
%   tab: table of morphological properties of soot using the benchmark...
%       ...correlations.
% ----------------------------------------------------------------------- %

%%% define constants %%%
dpp_100 = 17.8; % universal correlation prefactor
D_TEM = 0.35; % universal correlation exponent
k_a = 1.16; % screening prefactor
alpha_a = 1.1; % screening exponent
k_f = 1.3; % fractal prefactor
D_f = 1.78; % fractal dimension
rho_m = 1860; % estimate for continuum density of black carbon
k_m = 2.93; % mass-mobility prefactor
D_m = 2.48; % mass-mobility exponent
rho_0 = 1000; % normalized density for aerodynamic diameter
magg_100 = 2.93e-6; % aggregate mass at 100 nm mobility diameter
%%%%%%%%%%%%%%%%%%%%%%%%


da = exp(log(dpp/dpp_100) / D_TEM + log(100)); % projected area diameter (nm)

npp = k_a * (da ./ dpp).^alpha_a; % number of primaries (not rounded!)

Rg = exp(log(npp/k_f) / D_f + log(dpp/2)); % radius of gyration (nm)

magg1 = 1e18 * rho_m * pi * (npp .* (1e-9*dpp).^3) / 6; % mass obtained from geometry (fg)

magg2 = magg_100 * (dm / 100).^D_m; % mass from scaling law (fg)

dve = 1e9 * (6 * (1e-18*magg2) / (rho_m * pi)).^(1/3); % volume equivalent diameter (nm)

rho_eff = (6 * (1e18*magg2)) ./ (pi * (1e-9*dm).^3); % effective density (kg/m3)

x = (rho_m./rho_eff).^(1/3) .* (CUNNINGHAM(dve)./CUNNINGHAM(dm)); % shape factor

% trial and error for aerodynamic diameter
dae_i = 2 * rand(length(dm),1) .* dm; % initial guess
dae = dm;
ind = 0;
while any(abs((dae - dae_i) ./ dae_i) > 0.01) && (ind < 1e3)
    if ind~=0; dae_i = dae; end
    
    dae = dm.^((D_m - 1) / 2) .* sqrt((CUNNINGHAM(dm)./CUNNINGHAM(dae_i)) *...
        6 * k_m / (pi * rho_0));
    
    ind = ind + 1;
end

tab = table(dm, magg1, magg2, rho_eff, dve, x, dae, dpp,...
    npp, Rg); % make the table of properties

end

function cc = CUNNINGHAM(dp) % Cunningham correction factor

% constant coefficients
A1 = 1.257;
A2 = 0.4;
A3 = 0.55;

lambda = 6.8e-8; % ambient air mean free path

cc = 1 + (2 * lambda ./ dp) .* (A1 + A2 * exp(-A3*dp/lambda));

end
