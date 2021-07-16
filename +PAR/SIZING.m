function [dv, dg, dmax, dpp] = SIZING(pars)
% "SIZING" computes various equivalent sizes proposed for fractal...
%   ...aggregates.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   pars: Particle information structure/class
% ----------------------------------------------------------------------- %
%
% Output:
%   d_v: Equivalent volumetric (Smoluchowski's) diameter
%   d_g: Gyration diameter
%   d_max: % Maximum extent diameter
%   d_pp: % Mean + std of sizes of primary particles
% ----------------------------------------------------------------------- %

% Total number of (independent) particles
if isa(pars, 'AGG')
    n_par = size(pars, 1);
else
    n_par = size(pars.n, 1);
end

% Initializing the size arrays
dv = zeros(n_par,1);
dg = zeros(n_par,1);
dmax = zeros(n_par,1);
dpp = zeros(n_par,1);

for i = 1 : n_par
    
    if isa(pars, 'AGG')
        pp = AGG.COMPILEPP(pars(i)); % Compiling primary particle info
        n_pp = size(pp,1); % Number of primaries
        r_com = AGG.COM(pp); % Center of mass
        dmax(i) = AGG.TERRITORY(pp); % Maximum extent
        dpp(i,1) = mean(pars(i).pp.d); % Average diameter of primaries
        dpp(i,2) = std(pars(i).pp.d); % Sample standard deviation of...
            % ...primary particle size
        
    else
        pp = pars.pp{i};
        n_pp = size(pp,1);
        r_com = PAR.COM({pp}, n_pp);
        dmax(i) = PAR.TERRITORY({pp}, n_pp);
        dpp(i,1) = mean(pars.pp{i}(:,2));
        dpp(i,2) = std(pars.pp{i}(:,2));
    end
       
    dv(i) = nthroot(sum(pp(:,2) .^ 3), 3); % Volumetric diameter
    
    a2 = sum((repmat(r_com, n_pp, 1) - pp(:,3:5)).^2, 2); % Distance...
        % ...between center of primaries and aggregate center of mass
    radg2 = 0.6 .* (pp(:,2) ./ 2).^2; % Radius of gyration of...
        % ...primaries squared
    dg(i) = 2 .* sqrt(sum(pp(:,2).^3 .* (a2 + radg2)) ./ sum(pp(:,2).^3));
        % Gyration diameter
    
end

end

