function [pp, r_com] = TRANSLATE(pp, r_com, n_pp, dr)
% "TRANSLATE" moves the center of mass and primary particle locations of...
%     ...an aggregate via a translational motion.
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     r_com: Center of mass of (independent) particles
%     pp: primary particles info cell array
%     n_pp: Number of primaries within each aggregate
%     dr: Translation vector
% ----------------------------------------------------------------------- %

% Translating the primary particles
pp = cell2mat(pp);
pp(:,3:5) = pp(:,3:5) + repelem(dr, n_pp, 1);
pp = mat2cell(pp, n_pp);

% Moving the center of mass
r_com = r_com + dr;

end

