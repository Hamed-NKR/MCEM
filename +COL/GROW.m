 function [pars, timetable] = GROW(pars, timetable)
% "GROW" monitors collision of particles and updates the particle...
%     ...structure based on the post-aggregation data.
% ----------------------------------------------------------------------- %
% 
% Input/Output:
%     pars: Particle info structure/class
% ----------------------------------------------------------------------- %

% Total number of (independent) particles
if isa(pars, 'AGG')
    n_par = length(pars);
else
    n_par = size(pars.n, 1);
end

% Compiling/copying properties locally
dmax = cat(1, pars.dmax); % Maximum distance within the aggregates from...
    % ...their center of mass
r = cat(1, pars.r);

% timetable.preoverlap = [timetable.preoverlap; clock];
% Making particle pair indices
ind_pars = (1 : n_par)';
ind_pars = [repelem(ind_pars, n_par, 1), repmat(ind_pars, n_par, 1)];

% identifying repeating pairs
rmv1 = (1 : n_par)';
rmv1 = repelem((rmv1 - 1) .* n_par, 1 : n_par);
rmv2 = repmat((1 : n_par)', [1 n_par]);
rmv2 = triu(rmv2);
rmv2 = reshape(rmv2, n_par^2, 1);
rmv2(rmv2 == 0) = [];
rmv = rmv1 + rmv2;
ind_pars(rmv,:) = [];

% Generating the "OVR" inputs:
d_pairs = [repelem(dmax, n_par, 1), repmat(dmax, n_par, 1)];
    % size input
d_pairs(rmv,:) = [];
r_pairs = [repelem(r, n_par, 1), repmat(r, n_par, 1)]; % Pairs of...
    % ...coordinates
r_pairs(rmv,:) = [];

% Checking overlapping
ovrs = COL.OVR(r_pairs, d_pairs);
% timetable.postoverlap = [timetable.postoverlap; clock];

% Updating the location of overlapped particles
if ~ isempty(find(ovrs == 1, 1))
    
    ind_chk = ind_pars(ovrs == 1, :); % Indices of colliding pairs
    n_chk = size(ind_chk,1); % Number of colliding pairs
    
    for i = 1 : n_chk
        if ind_chk(i,1) ~= ind_chk(i,2)
            
%             timetable.preconnect = [timetable.preconnect; clock];
            % Sticking the two colliding particles
            if isa(pars, 'AGG')
                pp1 = AGG.COMPILEPP(pars(ind_chk(i,1)));
                pp2 = AGG.COMPILEPP(pars(ind_chk(i,2)));
                [pp1, pp2, colstat] = COL.CONNECT(pp1, pp2);
                pars(ind_chk(i,1)).r = pp1(:,3:5);
                pars(ind_chk(i,2)).r = pp2(:,3:5);
            else
                [pars.pp{ind_chk(i,1)}, pars.pp{ind_chk(i,2)}, colstat]...
                    = COL.CONNECT(pars.pp{ind_chk(i,1)},...
                    pars.pp{ind_chk(i,2)}); 
            end
%             timetable.postconnect = [timetable.postconnect; clock];
            
%             % Initializing collision status variable (1 --> collided,...
%                 % ...0 --> uncollided)
%             if ~exist('colstat', 'var'); colstat = 1; end
            
%             timetable.preunite = [timetable.preunite; clock];
            if colstat
                
                % Merging the particles info
                [pars, ind_new] = COL.UNITE(pars, [ind_chk(i,1),...
                    ind_chk(i,2)]);

                % Updating the general overlap checking indices after...
                    % ...each unification
                for j = 1 : n_chk
                    ind_chk(j,1) = ind_new(find(ind_new(:,1) ==...
                        ind_chk(j,1), 1), 2);
                    ind_chk(j,2) = ind_new(find(ind_new(:,1) ==...
                        ind_chk(j,2), 1), 2);
                end
            end
%             timetable.postunite = [timetable.postunite; clock];
            
        end
    end
end

end