function pp = INIT_MORPH(pp)
% "INIT_MORPH" assigns initial (random) morpholgy of the primaries...
%     ...within individual aggregates.
% 
% NOTE: Only suggested for initialization of aggregates with...
%     ...a few number of primaries.
% ----------------------------------------------------------------------- %
% 
% Input/Output:
% pp: An array of primary particle information (containing their...
%     ...indices, sizes, and locations)
% ----------------------------------------------------------------------- %

n_par = size(pp,1); % Total number of aggregates
n_pp = zeros(n_par, 1); % Number of primaries within the particles

for i = 1 : n_par
    
    n_pp(i) = size(pp{i},1);
    
    if n_pp(i) > 1
        
        sel_pp = randperm(n_pp(i)); % Random primary particle...
            % ...selection array
        
        for j = 2 : n_pp(i)
            
            chk = 1; % Overlap checking criterion   
            
            while ~ isempty(find(chk == 1, 1))
                
                jj1 = sel_pp(j); % New primary particle index
                jj0 = sel_pp(randperm(j-1,1)); % Finding a random index...
                    % ...corresponding to the target primary particle
                
                phi = 2 * pi * rand; % A random azimuthal orientation...
                    % ...for deposition of the next primary
                theta = pi * rand; % A random polar orientation for...
                    % ...deposition
                
                % Spherical coordinates for the new primary particle...
                    % ...center with respect to the target
                disp1 = [abs(pp{i}(jj1,2) + ...
                    pp{i}(jj0,2)) / 2; theta; phi];
                
                % Calculating the absolute center position for the...
                    % ...newly deposited primary particle
                pp{i}(jj1,3:5) = pp{i}(jj0,3:5) + ...
                    [disp1(1) * sin(disp1(2)) * cos(disp1(3)), ...
                    disp1(1) * sin(disp1(2)) * sin(disp1(3)), ...
                    disp1(1) * cos(disp1(2))];

                % Checking for overlapping with the previous primaries
                r_ovr = [repmat(pp{i}(jj1,3:5), j-1, 1), ...
                    pp{i}(sel_pp(1 : j-1),3:5)];
                d_ovr = [repmat(pp{i}(jj1,2), j-1, 1), ...
                    pp{i}(sel_pp(1 : j-1),2)];
                chk = COL.OVR(r_ovr, d_ovr);
                
            end
            
            % Post-deposition intrinsic rotation
            angs = 2 * pi * rand(3,1); % A set of 3 random Euler angles...
                % ...(yaw, pitch, and roll)
            pp(i) = PAR.ROTATE(pp(i), n_pp(i), angs);
            
        end
        
    end
    
end

end
