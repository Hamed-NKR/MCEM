
% AGG  A class containing the aggregate properties and functions.
%  
%  05-2021

classdef AGG
    
    properties
        pp = struct()   % primary particle information
        n = [];         % number of primary particles
        
        dmax = []       % maximum extent of diameter
        r = []          % center of mass
        v = []          % aggregate velocity
        
        m = []          % ~ mass
        rho = []        % ~ effective density
        delt = []       % ~ motion time-step
        tau = []        % ~ relaxation time
        f = []          % ~ friction factor
        diff = []       % ~ diffusivity
        lambda = []     % ~ diffusive mean free path
        kn = []         % Knudsen number (both kinetic and diffusive)
        nnl = []        % ~ nearest neighbor list
    end
    
    methods
        %== AGG ==========================================================%
        function obj = AGG(pp)
         % Function to initialize the aggregate object.
            
            if nargin == 0; return; end
            
            % Assign aggregate velocity.
            if ~exist('v', 'var'); v = []; end
            if isempty(v); v = [0, 0, 0]; end  % no velocity be default
            obj.v = v;
            
            % Assign primary particle information.
            obj.pp = pp;
            obj.n = size(pp.r, 1);
            
            obj.r = obj.COM(pp);
            obj.dmax = obj.TERRITORY(pp, obj.r);
            
        end
        
        
        %== TRANSLATE ====================================================%
        function objs = TRANSLATE(objs, dr)
            
            for ii=1:length(objs)
                objs(ii).pp.r = objs(ii).pp.r + dr(ii, :);
                
                objs(ii).r = objs.COM(objs(ii).pp);
            end
            
            % TERRITORY (max. outer diameter) does not change.
        end
        
        
        %== RENDER =======================================================%
        function h = RENDER(obj, idx, cm)
            
            % Select only idx (first if not given) aggregate.
            if ~exist('idx', 'var'); idx = []; end
            if isempty(idx); idx = 1; end
            obj = obj(idx);
            
            % Set default colormap as "summer".
            if ~exist('cm', 'var'); cm = []; end
            if isempty(cm); cm = summer; end
            
            pp0 = obj.pp;  % shorten subsequent references to pp
            n_pp = length(pp0.dp);
            
            disp('Rendering:');
            UTILS.TEXTBAR([0, n_pp]);
            
            clf; axis equal; hold on;
            colormap(cm);
            
            % Plot spheres.
            [X,Y,Z] = sphere(60);
            for ii=1:n_pp
                h = surf(X .* pp0.dp(ii) ./ 2 + pp0.r(ii,1), ...
                    Y .* pp0.dp(ii) ./ 2 + pp0.r(ii,2), ...
                    Z .* pp0.dp(ii) ./ 2 + pp0.r(ii,3));
                lightangle(-45,30)
                h.FaceLighting = 'gouraud';
                h.AmbientStrength = 0.8;
                h.DiffuseStrength = 0.2;
                h.SpecularStrength = 0.05;
                h.SpecularExponent = 2;
                h.BackFaceLighting = 'lit';
                
                UTILS.TEXTBAR([ii, n_pp]);
            end
            disp(' ');
            
            % Format plot.
            disp('Formatting plot ...');
            camlight('right');
            shading interp;
            view([-37, 20]);
            h = gca;
            axis('off');
            
            disp('DONE.');
            disp(' ');
            
            figure(gcf);
            hold off;
            
            if nargout==0
                clear h; 
            end
        end
        
    end
    
    
    
    methods(Static)
        
        %== COM ======================================================%
        function com = COM(pp)
        % Computed the center of mass of the aggregate.
            
            m = pp.dp .^ 3;
            
            com = pp.r .* m;
            com = sum(com, 1) ./ sum(m);
        end
        
        
        %== GET_MAX ======================================================%
        function dmax = TERRITORY(pp, com)
        % Get the maximum extent of the aggregate.
            
            % Most distant primary ppticle.
            [dmax, imax] = ...
                max(sum(sqrt((pp.r - com) .^ 2), 2));
            
            dmax = dmax + pp.dp(imax) ./ 2;  % add the primary ppt. diameter
            
        end
        
        
        %== CHECK_OVERLAP ===========================================s=====%
        function f_overlap = OVR(pp)
        % Output is a logical indicating if a ppticle is overlapping another.
            
            % Loop through ppticles to check overlap.
            
        end
        
        
        %== COMPILEPP ====================================================%
        %   Compile primary particle information across multiple
        %   aggregates.
        function pp = COMPILEPP(objs)
            
            pp = [];
            for ii=1:length(objs)
            	pp = [pp; ...
                      [objs(ii).pp.id, ...
                       objs(ii).pp.dp, ...
                       objs(ii).pp.r]];
            end
        end
        
        %== COMPILE ======================================================%
        %   Compile another property over multiple aggregates.
        function prop = COMPILE(objs, propname)
            
            prop = [];
            for ii=1:length(objs)
            	prop = [prop; objs(ii).(propname)];
            end
        end
    end
end

