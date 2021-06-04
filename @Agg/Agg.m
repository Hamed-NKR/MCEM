
% AGG  A class containing the aggregate properties and functions.
%  
%  05-2021

classdef Agg
    
    properties
        pp = struct()   % primary particle information
        dmax = []       % maximum extent of diameter
        com = []        % center of mass
    end
    
    methods
        %== AGG ==========================================================%
        function obj = Agg(pp)
         % Function to initialize the aggregate object.
            
            obj.pp = pp;
            
            obj.com = obj.GET_COM(pp);
            obj.dmax = obj.GET_MAX(pp, obj.com);
            
        end
        
        
        %== TRANSLATE ====================================================%
        function obj = TRANSLATE(obj, dr)
            obj.pp.r = obj.pp.r + dr;
            
            obj.com = obj.GET_COM(obj.pp);
            % DMAX does not change.
        end
        
        
        function h = RENDER(obj)
            pp0 = obj.pp;  % shorten subsequent references to pp
            n_pp = length(pp0.dp);
            
            disp('Rendering:');
            UTILS.TEXTBAR([0, n_pp]);
            
            clf; axis equal; hold on;
            colormap(summer);
            
            % Plot spheres.
            [X,Y,Z] = sphere(60);
            for ii=1:n_pp
                h = surf(X .* pp0.dp(ii) + pp0.r(ii,1), ...
                    Y .* pp0.dp(ii) + pp0.r(ii,2), ...
                    Z .* pp0.dp(ii) + pp0.r(ii,3));
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
        end
        
    end
    
    
    
    methods(Static)
        
        %== GET_COM ======================================================%
        function com = GET_COM(pp)
        % Computed the center of mass of the aggregate.
            
            m = pp.dp .^ 3;
            
            com = pp.r .* m;
            com = sum(com, 1) ./ sum(m);
        end
        
        
        %== GET_MAX ======================================================%
        function dmax = GET_MAX(pp, com)
        % Get the maximum extent of the aggregate.
            
            % Most distant primary particle.
            [dmax, imax] = ...
                max(sum(sqrt((pp.r - com) .^ 2), 2));
            
            dmax = dmax + pp.dp(imax) ./ 2;  % add the primary part. diameter
            
        end
        
        
        %== CHECK_OVERLAP ===========================================s=====%
        function f_overlap = CHECK_OVERLAP(pp)
        % Output is a logical indicating if a particle is overlapping another.
            
            % Loop through particles to check overlap.
            
        end
    end
end

