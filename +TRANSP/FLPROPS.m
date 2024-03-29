function [mu_f, lambda_f] = FLPROPS(fl, params_const, opts)
% "KINETIC" calculates the fluid slip-related properties based on the...
%     ...kinetic theory of gases.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     fl: Fluid information structure
%     params_const: Problem's table of constant physical properties
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     mu_f: Fluid viscosity
%     lambda_f: Fluid mean free path
% ----------------------------------------------------------------------- %

% M_air = params_const.Value(2); % Air molar mass (kg/mol)
% Na = params_const.Value(4); % Avogadro constant (mol^-1)
% Ru = params_const.Value(5); % Universal gas constant (j/mol.k)
% 
% d_kin = (0.79 * 3.64 + 0.21 * 3.46) * 1e-10; % Air kinetic diameter...
% % ...(an average estimate based on values for N2 and O2)
% sig_col = pi * (d_kin^2); % Collision area (m2)
% mu_f = 2 * sqrt(M_air*Ru*(fl.temp)) / (3*sqrt(pi)*Na*sig_col); % Viscosity
% lambda_f = mu_f / (0.499*(fl.p)*sqrt(8*M_air/(pi*Ru*(fl.temp)))); % Mean...
%     % ...free path 

% initialize ambient variable
if ~exist('opts', 'var') 
    opts = struct();
end
if (~isfield(opts, 'amb')) || isempty(opts.amb)
    opts.amb = 'flame'; % assume flame conditions as default
end

if strcmp(opts.amb, 'room') % @ 20 c & 1 atm: room conditions; typical of post-flame agglomeration
    mu_f = 1.813e-5; % Viscosity (Pa.s)
    lambda_f = 6.8e-8; % Mean free path (m)
elseif strcmp(opts.amb, 'flame') % @ 1700 k & 1 atm: typical of soot formation within flame
    mu_f = 5.647e-5;
    lambda_f = 4.98e-7;
end

end

