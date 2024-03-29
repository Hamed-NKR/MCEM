function [pars, timetable] = GROW(pars, opts, timetable)
% "GROW" monitors collision of particles and updates the particle...
%     ...structure based on the post-aggregation data.
% ----------------------------------------------------------------------- %
% 
% Input/Output:
%     pars: Particle info structure/class
% ----------------------------------------------------------------------- %

% initialize option input if not given
if ~exist('opts', 'var')
    opts = struct();
end

if (~isfield(opts, 'indupdate')) || isempty(opts.indupdate)
    opts.indupdate = 'on'; % default to update agg indices after collision
end

opts_indupdate = opts.indupdate;

% set default for collision mechanism
if (~isfield(opts, 'col')) || isempty(opts.col)
    opts.col = 'agg'; % default to be aggregation/agglomeration, not coalescense
end

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

ind_pairs = nchoosek(1 : n_par, 2); % Making particle pair indices
d_pairs = [dmax(ind_pairs(:,1)), dmax(ind_pairs(:,2))]; % Pair of sizes
r_pairs = [r(ind_pairs(:,1),:), r(ind_pairs(:,2),:)]; % Pair of coordinates

% Checking overlapping
ovrs = COL.OVR(r_pairs, d_pairs);

% timetable.postoverlap = [timetable.postoverlap; clock];

% Updating the location of overlapped particles
if any(ovrs == 1)
    
    ind_chk = ind_pairs(ovrs == 1, :); % Indices of colliding pairs
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
                    ind_chk(i,2)], opts.col);

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
    
    if ismember(opts_indupdate, {'ON', 'On', 'on'})
        for k = 1 : length(pars.n)
            pars.pp{k}(:,6) = k;
        end
    end
end


end