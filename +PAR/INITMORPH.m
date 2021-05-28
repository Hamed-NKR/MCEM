function [par_pp, par_d] = INITMORPH(par_pp, par_n)
% "INITMORPH" assigns initial (random) morpholgy of the primaries within...
    % individual aggregates.

% NOTE: Only suggested for initialization of aggregates with...
    % ...a few number of primaries.

% Inputs/Outputs:
% par_pp: An array of primary particle information (containing their...
    % ...indices, sizes, and locations)
% par_n: Size distribution of primaries
% par_r: Spatial location of the (independent) particle centers

n_par = size(par_pp,1); % Total number of aggregates

for i = 1 : n_par
    
    n_pp = size(par_pp{i},1); % Number of primaries within the particles
    
    if n_pp > 1
        
        sel_pp = randperm(n_pp); % Random primary particle...
            % ...selection array
        
        for j = 2 : n_pp
            
            chk = 1; % Overlap checking criterion   
            
            while chk
                
                jj1 = sel_pp(j); % New primary particle index
                jj0 = sel_pp(randperm(j-1,1)); % Finding a random index...
                    % ...corresponding to the target primary particle
                
                phi = 2 * pi * rand; % A random azimuthal orientation...
                    % ...for deposition of the next primary
                theta = pi * rand; % A random polar orientation for...
                    % ...deposition
                
                % Spherical coordinates for the new primary particle...
                    % ...center with respect to the target
                disp1 = [abs(par_pp{i}(jj1,2) + ...
                    par_pp{i}(jj0,2)) / 2; theta; phi];
                
                % Calculating the absolute center position for the...
                    % ...newly deposited primary particle
                par_pp{i}(jj1,3:5) = par_pp{i}(jj0,3:5) + ...
                    [disp1(1) * sin(disp1(2)) * cos(disp1(3)), ...
                    disp1(1) * sin(disp1(2)) * sin(disp1(3)), ...
                    disp1(1) * cos(disp1(2))];

                % Checking for overlapping with the previous primaries
                r_ovr = [repmat(par_pp{i}(jj1,3:5), j-1, 1), ...
                    par_pp{i}(sel_pp(1 : j-1),3:5)];
                d_ovr = [repmat(par_pp{i}(jj1,2), j-1, 1), ...
                    par_pp{i}(sel_pp(1 : j-1),2)];
                chk = COL.OVR(r_ovr, d_ovr);
                
            end
            
            % Post-deposition intrinsic rotation
            angs = 2 * pi * rand(3,1); % A set of 3 Euler angles...
                % ...(yaw, pitch, and roll)
            yaw = [cos(angs(1)), -sin(angs(1)), 0;...
                sin(angs(1)), cos(angs(1)), 0; 0, 0, 1];
                % transformation matrix for yaw rotation
            pitch = [cos(angs(2)), 0, sin(angs(2));...
                 0, 1, 0; -sin(angs(2)), 0, cos(angs(2))];
                % transformation matrix for pitch rotation
            roll = [1, 0, 0; 0, cos(angs(3)), -sin(angs(3));...
                0, sin(angs(3)), cos(angs(3))];
                % transformation matrix for roll rotation
            rot = yaw * pitch * roll; % The net rotation matrix
            % Obtaining the post rotation coordinates of primaries
            par_pp{i}(:,3:5) = (rot * (par_pp{i}(:,3:5))')';
            
        end
        
    end
    
end

[~, par_d] = COL.EQUIV(par_pp, par_n); % Assigning radius of gyration as...
    % ...the equivalent size of the aggregates

end
