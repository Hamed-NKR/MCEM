function [par_pp, par_r] = TRANSLATE(par_pp, par_r, par_n, dr)
% "TRANSLATE" moves the center of mass and primary particle locations of...
%     ...an aggregate via a translational motion.
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     par_r: Center of mass of (independent) particles
%     par_pp: primary particles info cell array
%     par_n: Number of primaries within each aggregate
%     dr: Translation vector
% ----------------------------------------------------------------------- %

% Translating the primary particles
pp = cell2mat(par_pp);
pp(:,3:5) = pp(:,3:5) + repelem(dr, par_n, 1);
par_pp = mat2cell(pp, par_n);

% Moving the center of mass
par_r = par_r + dr;

end

