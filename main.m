%% General description

% The present script:
% 1. reads the aggregation input parameters,
% 2. randomly initializes the monomers' locations and velocities,
% 3. moves the particles based on brownian motions and drag,
% 4. imposes periodic conditions on the boundaries,
% 5. checks for the particles' collisions over time,
% 6. merge the collided particles as larger clusters,
% 7. maintians the volume fraction by enlarging the main domain upon
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

% Defining the domain structure (for possible later uses)
dom = struct('size',[]);
dom.size = dom_size;
% Input is the domain size array.

% Loading the fluid properties
fl = struct('temp',temp_f,'v',v_f,'p',p_f,'mu',[],'lambda',[]);
% fl is the fluid information structure for particle-fluid interactions.
[fl.mu, fl.lambda] = CPL.KINETIC(fl.temp,fl.p); % fluid viscosity and...
% ...mean free path

% Declaring the primary particles structure
pp = struct('ind',[1:n_pp]','ind_agg',zeros(n_pp,1),'d',[],...
    'r',[],'v',[],'nn',[]);
% Inputs are primary particle index, corresponding aggregate index,...
% position and velocity of the primaries, and their nearest neighbor list.
% pp rows correspond to different primaries.

% Assigning the primary particle diameters
pp.d = PP.INIT.DIAM(n_pp,d_pp);

% Assigning the primary particle initial locations
pp.r = PP.INIT.LOC(dom_size,pp.d);

% Assigning the primary particle initial velocities
pp.v = PP.INIT.VEL(pp.d,fl.temp);

disp("The computational domain is successfully initialized...")

fig_pp_init = figure(1);
VIS.PLOTPP(dom_size,pp,1)

fig_pp_nn = figure(2);
VIS.PLOTNN(dom_size,pp,randperm(n_pp,1))

%% Solving equation of motion for the particles

k_max = 300; % Marching index limit
time = zeros(k_max,1);

fig_pp_anim = figure(3);
disp('Simulating:');
UTILS.TEXTBAR([0, k_max]);  % Initializing textbar
UTILS.TEXTBAR([1, k_max]);  % Indicating start of marching

for k = 2 : k_max
    
    [pp.r, pp.v, delt] = MOV.MARCH(pp,fl); % Solving equation of motion
    pp.r = MOV.PBC(dom_size,pp.r); % Applying periodic boundary conditions
    
    t_plot = 1;
    if mod(k-1,t_plot) == 0
        VIS.PLOTPP(dom_size,pp,0) % Plotting every t_plot time steps
        drawnow; % Drawing the plot at the desired time steps
    end
    
    time(k) = time(k-1) + delt; % Updating time (this needs to be inside...
    % the loop for later uses)
    
    UTILS.TEXTBAR([k, k_max]);  % Updating textbar
    
end
