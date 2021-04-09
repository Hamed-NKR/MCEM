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
    'r',[],'v',[],'delt',[],'tau',[],'diff',[],'lambda',[],'nn',[]);
% Inputs are primary particle index, corresponding aggregate index,...
% position and velocity of the primaries, their diffusive properties,...
% ,and nearest neighbor list.
% Rows of pp elements correspond to different primaries information.

% Assigning the primary particle diameters
pp.d = PP.INIT.DIAM(n_pp,d_pp);

% Assigning the primary particle initial locations
pp.r = PP.INIT.LOC(dom_size,pp.d);

% Assigning the primary particle initial velocities
pp.v = PP.INIT.VEL(pp.d,fl.temp);

disp("The computational domain was successfully initialized...")

fig_pp_init = figure(1);
VIS.PLOTPP(dom_size,pp,1);

fig_pp_nn = figure(2);
VIS.PLOTNN(dom_size,pp,randperm(n_pp,1),10);

%% Solving equation of motion for the particles

k_max = 1000;  % Marching index limit
time = zeros(k_max,1);
fig_pp_anim = figure(3);
t_plt = 1;  % Defining a plotting timeframe criterion

prompt = 'Do you want the animation to be saved? Y/N: ';
str = input(prompt,'s');
if (str == 'Y') || (str == 'y')
    video_pp = VideoWriter('outputs\Animation_DLCA.avi');  % Initializing video
    video_pp.FrameRate = 5;  % Setting frame rate
    open(video_pp);  % Opening video file
end

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]);  % Initializing textbar
UTILS.TEXTBAR([1, k_max]);  % Indicating start of marching

for k = 2 : k_max
    
    [pp, delt] = MOV.MARCH(pp, fl); % Solving equation of motion
    pp.r = MOV.PBC(dom_size,pp.r); % Applying periodic boundary conditions
    
    if mod(k-2,t_plt) == 0
        VIS.PLOTPP(dom_size, pp, 0); % Plotting every t_plt time steps
        drawnow; % Drawing the plot at the desired time steps
        pause(0.1); % Slowing down the animation speed
        if (str == 'Y') || (str == 'y')
            frame_now = getframe(fig_pp_anim);  %  Capturing current frame
            writeVideo(video_pp, frame_now);  %  Saving the video
        end
    end
    
    time(k) = time(k-1) + delt; % Updating time (this needs to be inside...
    % the loop for later uses)
    
    UTILS.TEXTBAR([k, k_max]);  % Updating textbar
    
end

if (str == 'Y') || (str == 'y')
    close(video_pp);  % Closing the video file
end


