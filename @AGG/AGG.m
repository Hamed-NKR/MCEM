classdef AGG
% "AGG" is a class containing the aggregate properties and functions.
%
% Original author: Timothy Sipkens, 05-2021
% Revised by: Hamed Nikookar, 07-2021
% ----------------------------------------------------------------------- %
    
    properties
        
        pp = struct('r', [], 'd', [], 'id', [])   % Primary particle...
            % ...properties structure
        n = []          % Number of primary particles within the aggregates
        
        dv = []         % Aggregates equivalent volumetric diameter
        dg = []         % ~ gyration diameter
        dm = []         % ~ mobility diameter
        da = []         % ~ aerodynamic diameter
        dpp = []        % ~ mean + std of diameter of primaries
        dmax = []       % ~ maximum extent of diameter
        
        r = []          % ~ center of mass
        v = []          % ~ velocity
        
        m = []          % ~ mass
        rho = []        % ~ effective density
        delt = []       % ~ motion time-step
        tau = []        % ~ relaxation time
        f = []          % ~ friction factor
        diff = []       % ~ diffusivity
        lambda = []     % ~ diffusive mean free path
        kn_kin = []     % ~ kinetic Knudsen number
        kn_diff = []    % ~ diffusive Knudsen number
        nnl = []        % ~ nearest neighbor list
        
    end
% ----------------------------------------------------------------------- %
    
    methods
        
        % === AGG ======================================================= %
        function obj = AGG(pp)
        % "AGG" initializes the aggregate objects.
        % --------------------------------------------------------------- %
        %
        % pp: Primary particle structure
        % obj: Aggregate object
        % --------------------------------------------------------------- %
            
            if nargin == 0; return; end
            
            % Assigning the aggregate properties
            obj.pp = pp;
            obj.n = size(pp.id, 1);
            obj.r = obj.COM(pp);
            if isempty(obj.v); obj.v = [0, 0, 0]; end  % Zero velocity...
                % ...be default
        end
        
        % === TRANSLATE ================================================= %
        function objs = TRANSLATE(objs, dr)
        % "TRANSLATE" moves the aggregates in certain directions.
        % --------------------------------------------------------------- %
        %
        % objs: Aggregate objects
        % dr: translation vectors
        % --------------------------------------------------------------- %
        
            for i = 1 : length(objs)
                objs(i).pp.r = objs(i).pp.r + dr(i,:); % Moving the...
                    % ...primaries
                objs(i).r = objs(i).r + dr(i,:); % ~ the center of mass
            end
            
        end
        
        % === ROTATE ==================================================== %
        function objs = ROTATE(objs, dtheta)
        % "ROTATE" turns the aggregates in certain intrinsic angles.
        % --------------------------------------------------------------- %
        %
        % objs: Aggregate objects
        % dtheta: Intrinsic rotation angles
        % --------------------------------------------------------------- %
            
            n_objs = length(objs); % Number of aggregates
            
            for i = 1 : n_objs
                
                % Transformation matrix for yaw rotation
                yaw = [cos(dtheta(i,1)), -sin(dtheta(i,1)), 0;...
                    sin(dtheta(i,1)), cos(dtheta(i,1)), 0; 0, 0, 1];
                    
                % ~ pitch
                pitch = [cos(dtheta(i,2)), 0, sin(dtheta(i,2));...
                     0, 1, 0; -sin(dtheta(i,2)), 0, cos(dtheta(i,2))];
                
                % ~ roll
                roll = [1, 0, 0; 0, cos(dtheta(i,3)), -sin(dtheta(i,3));...
                    0, sin(dtheta(i,3)), cos(dtheta(i,3))];

                % The net rotation matrix
                rot = yaw * pitch * roll;
                
                % Rotating pimaries around the center of mass
                r_com = objs(i).COM(objs(i).pp);
                objs(i).pp.r = (rot * (objs(i).pp.r - r_com)')' + r_com;
                
            end
        end
        
    end
% ----------------------------------------------------------------------- %
    
    methods(Static)
        
        % === COM ======================================================= %
        function r_com = COM(pp)
        % "COM" computes the center of mass of an ensemble of primaries.
        % --------------------------------------------------------------- %
        %
        % pp: Primary particle structure
        % r_com: Center of mass coordinates
        % --------------------------------------------------------------- %
            
            c_m = pp.d .^ 3; % Volumetric coefficients
            r_com = pp.r .* c_m;
            r_com = sum(r_com, 1) ./ sum(c_m);
            
        end
        
        % === TERRITORY ================================================= %
        function dmax = TERRITORY(pp)
        % "TERRIROY" gets the maximum extent of the aggregate.
        % --------------------------------------------------------------- %
        %
        % pp: Primary particle structure
        % dmax: Maximum distance of aggregate domain from its center of...
        %   ...mass
        % --------------------------------------------------------------- %
            
            r_com = AGG.COM(pp); % Center of mass coordinates
            dmax = max(sum(sqrt((pp.r - r_com) .^ 2) + pp.d ./ 2, 2));
            
        end
        
        % === COMPILEPP ================================================= %
        % "COMPILEPP" compiles a primary particle information set across...
        %   ...multiple aggregates.
        %
        % Note: This will be useful for some functions that aim to be...
        %   ...able to run in both "structure" and "class" modes.
        % --------------------------------------------------------------- %
        %
        % objs: Aggregate objects
        % pp_glob: Concatinated primary particle data
        % --------------------------------------------------------------- %
        
        function pp_glob = COMPILEPP(objs)
            
            n = cat(1, objs.n); % Compiling number of primaries data
            n_pp = sum(n); % Total number of primaries
            pp_glob = zeros(n_pp, 5); % Initializing the global array
            ii = 1; % Data calling index
            
            for i = 1 : length(objs)
                % Adding the i^th primary's data
            	pp_glob(ii : ii - 1 + n(i), :) =...
                    [objs(i).pp.id, objs(i).pp.d, objs(i).pp.r];
                ii = ii + n(i); % Updating the number of row to be...
                    % ...called
            end
            
        end
        
        % === COMPILE =================================================== %
        % "COMPILEPROP" compiles a general property over multiple...
        %   ...aggregates.
        %
        % Note: The purpose is again to make different functions...
        %   ...compatible with both "structure" and "class" modes.
        % --------------------------------------------------------------- %
        %
        % objs: Aggregate objects
        % propname: Property to be concatinated
        % prop: The concatinated outcome
        % --------------------------------------------------------------- %
        
        function prop = COMPILEPROP(objs, propname)
            
            % Initializing the data storage array
            prop = zeros(length(objs), size(objs(1).(propname), 2));
            
            for i = 1 : length(objs)
            	prop(i,:) = objs(i).(propname); % Adding the i^th...
                    % ...aggregate data
            end
            
        end
        
    end
    
end

