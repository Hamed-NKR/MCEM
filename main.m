%% Description

% This script:
% 1. reads the aggregation input parameters,
% 2. randomly initializes the monomers' locations and velocities,
% 3. moves the particles based on brownian motions and drag,
% 4. imposes periodic conditions on the boundaries,
% 5. checks for the particles' collisions over time,
% 6. merge the collided particles as larger clusters,
% 7. maintians the volume fraction by emlargin the main domain upon
% clusteraions,
% 8. exports and plots the aggregation data.

clc
clear
close all

%% Importing the user-defined parameters

% Extracting data from the input file
file_id = fopen('inputs\MCEM_Params.txt','r');
params = textscan(file_id,'%f%s%s%s','Headerlines',2,...
    'Delimiter','\t','MultipleDelimsAsOne',1);
fclose(file_id);

params_val = cell2mat(params(1));
dom_size = params_val(1:3); % Domain dimensions (m)
n_pp = params_val(4); % Number of primaries
d_pp = params_val(5:6); % Mean size and standard deviation of primaries (m)
temp_f = params_val(7); % Flow temperature (k)
v_f = params_val(8:10); % Flow velocity vector (m/s)
p_f = params_val(11); % Flow pressure (pa)

%% Initializing the computational domain parameters

% Defining the domain structure
dom = struct('size',[]);
dom.size = dom_size;
% Input is the domain size array.

% Loading the fluid properties
fl = struct('temp',temp_f,'v',v_f,'p',p_f);
% fl is the fluid information structure for particle-fluid interactions.

% Declaring the primary particles structure
pp = struct('ind',[1:n_pp]','ind_agg',zeros(n_pp,1),'d',[],...
    'r',[],'v',[]);
% Inputs are primary particle index, corresponding aggregate index, and...
% position and velocity of the primaries.
% pp rows correspond to different primaries.

% Assigning the primary particle diameters
pp.d = PP.INIT.DIAM(n_pp,d_pp);

% Assigning the primary particle initial locations
pp.r = PP.INIT.LOC(dom_size,pp.d);

% Assigning the primary particle initial velocities
pp.v = PP.INIT.VEL(pp.d,fl.temp);

disp("The computational domain is successfully initialized...")

fig_pp_init = figure(1);
VIS.PLOTPP(dom.size,pp,1)

%% Solving equation of motion for the particles

k_max = 10000; % Marching index limit
time = zeros(k_max,1);
[delt_pp, tau_pp, mu_f, lambda_f] = MOV.PROPS(pp.d,fl); 
% Marching time step, primaries characteristic time, fluid viscosity and...
% mean free path

fig_pp_anim = figure(2);
disp('Simulating:');
UTILS.TEXTBAR([0, k_max]);  % initializing textbar
UTILS.TEXTBAR([1, k_max]);  % indicating start of marching

for k = 2 : k_max    
    % TO DO: There are probably too many inputs here. Can this be simplified?
    [pp.r, pp.v] = MOV.MARCH(pp,delt_pp,fl);
    % CHANGE: Assign r and v directly. The multi-step is slower and
    % wastes memory. 
    
    % CHANGE: You should plot after you update the results, or the last
    % step is not plotted. Especially given the VIS.STCPLOT at the
    % beginning plots the original points.
    % CHANGE: Plot every t_plot time steps. Change to t_plot = 1 to plot
    % every step (much slower).
    t_plot = 10;
    if mod(k-1, t_plot)==0
        VIS.PLOTPP(dom.size,pp,0)
        drawnow;  % draw the plot each time step
    end
    
    % CHANGE: I added this textbar utility to indicate progress in command
    % line.
    UTILS.TEXTBAR([k, k_max]);  % update textbar
    
    % TO DO: The particles seem to disappear over time. Likely a bug in
    % applying periodic boundary conditions. 
end
