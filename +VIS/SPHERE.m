function SPHERE(x_c,y_c,z_c,r)
% This function plots a sphere using its center and radius.

% Generating a unit sphere
[x,y,z] = sphere;

% Scaling the radius
x = r * x;
y = r * y;
z = r * z;

surf(x+x_c,y+y_c,z+z_c)

% % Approach 2:
%
% % Definition of transformation angles
% theta = 0 : pi/50 : 2*pi;
% phi = 0 : pi/50 : 2*pi;
% n_theta = size(theta,2);
% n_phi = size(phi,2);
% 
% % Spherical to Cartesian coordinate transformation
% x = (r * (cos(theta))) * sin(phi)' + x_c * ones(n_theta,n_phi);
% y = (r * (sin(theta))) * sin(phi)' + y_c * ones(n_theta,n_phi);
% z = (r * ones(n_theta)) * cos(phi)' + z_c * ones(n_theta,n_phi);
% 
% Location array organization
% x = x(:);
% y = y(:);
% z = z(:);
%
% surf(x,y,z)

end

