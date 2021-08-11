function pars = SIZING(pars)
% "SIZING" computes various equivalent sizes proposed for fractal...
%   ...aggregates.
% ----------------------------------------------------------------------- %
%
% Input/output:
%   pars: Particle information structure/class
% ----------------------------------------------------------------------- %
    
if isa(pars, 'AGG')
    for i = 1 : length(pars)
        pars(i).dv = AGG.EQUIV(pars(i).pp); % Equivalent volumetric...
            % ...diameter
        pars(i).dg = AGG.GYRATION(pars(i).pp); % Gyration diameter
        pars(i).dmax = AGG.TERRITORY(pars(i).pp); % Maximum extent diameter
        pars(i).dpp = AGG.MEANPP(pars(i).pp); % Average primary particle...
            % ...diameter
    end

else
    pars.dv = PAR.EQUIV(pars.pp);
    pars.dg = PAR.GYRATION(pars.pp, pars.n);
    pars.dmax = PAR.TERRITORY(pars.pp, pars.n);
    pars.dpp = PAR.MEANPP(pars.pp);
end

end

