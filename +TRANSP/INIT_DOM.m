function [pars, fl] = INIT_DOM(params_ud, params_const)
% "INIT_DOM" sets different fields of the domain data structures.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   params_ud: Table of constant parametrs
%   params_const: ~ user-defined parameters
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     pars: Data structure of particles
%     fl: ~ background fluid
% ----------------------------------------------------------------------- %

% Declaring the fluid structure and loading the corresponding imported...
    % ...properties
fl = struct('size', params_ud.Value(2:4), 'temp', params_ud.Value(11),...
    'v', params_ud.Value(12:14), 'p', params_ud.Value(15), 'mu', [],...
    'lambda', []);
% "fl" is a structure containing the main physical properties of the...
    % ...background fluid. It is used to consider and study the...
    % ...particle-fluid interactions.
% Fields:
    % size: Computational domain size
    % temp: Fluid temperature
    % v: ~ velocity
    % p: ~ pressure
    % mu: ~ viscosity
    % lambda: ~ mean free path

[fl.mu, fl.lambda] = TRANSP.FLPROPS(fl,params_const); % Calculating the...
    % ...fluid viscosity and mean free path

% Declaring the particle structure
pars = struct('pp', [], 'n', [], 'dv', [], 'dg', [], 'dm', [], 'da', [],...
    'dmax', [], 'dpp', [], 'dpp_g', [], 'r', [], 'v', [], 'm', [],...
    'rho', [], 'delt', [], 'tau', [], 'f', [], 'diff', [], 'lambda', [],...
    'kn_kin', [], 'kn_diff', [], 'nnl', [], 'n_hyb', []);
% "pars" is a structure conating the main physical properties of...
    % ...independent particles (whether being single monomers or...
    % ...aggregates) as well as the properties of their constituent...
    % ...primary particles (if applicable).
% Fields:
    % pp: Primary particle characteristics (an N*5 matrix containing...
        % ...their index, size, and location)
    % n: Number of primaries within the particles
    % dv: Particles equivalent volumetric diameter
    % dg: ~ gyration diameter
    % dm: ~ mobility diameter
    % da: ~ aerodynamic diameter
    % dmax: ~ maximum extent diameter
    % dpp: ~ mean + std of diameter of primaries
    % dpp_g: ~ geometric mean + std of ~
    % r: ~ spatial location
    % v: ~ velocity
    % m: ~ mass
    % rho: ~ effective density
    % delt: ~ motion time-step
    % tau: ~ relaxation time
    % f: ~ friction factor
    % diff: ~ diffusivity
    % lambda: ~ diffusive mean free path
    % kn_kin: Kinetic Knudsen number
    % kn_diff: Diffusive Knudsen number
    % nnl: ~ nearest neighbor list
    % n_hyb: number of monodispersity regions (in polydisperse hybrid aggs)
% NOTE: The rows of each field correspond to characteristics of each...
    % ...independent particle.

end

