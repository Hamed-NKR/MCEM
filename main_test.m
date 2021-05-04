% This script is written to test and monitor the performance and outputs of 
% the program. Some parts, as a result, are not the most efficient in terms
% of computational cost.

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

% Declaring the particle structure
par = struct('pp',[],'d',[],'r',[],'v',[],...
'delt',[],'tau',[],'diff',[],'lambda',[],'nnl',[]);
% Inputs are list of primaries characteristics...
% ...(index, size and local coordinates), particles' equivalent position...
% ...and velocity, their diffusive properties, and nearest neighbor list.
% Element rows correspond to different aggregates info.

% Assigning the primary particle diameters
par.d = PAR.INITDIAM(n_pp,d_pp);

% Assigning the primary particle initial locations
par.r = PAR.INITLOC(dom_size,par.d);

% Assigning the primary particle initial velocities
par.v = PAR.INITVEL(par.d,fl.temp);

disp("The computational domain was successfully initialized...")

fig_init = figure(1);
UTILS.PLOTPAR(dom_size,par,1);

fig_nn = figure(2);
UTILS.PLOTNN(dom_size,par,randperm(n_pp,1),10);

%% Solving equation of motion for the particles

k_max = 50;  % Marching index limit
time = zeros(k_max,1);
fig_anim = figure(3);
t_plt = 1;  % Defining a plotting timeframe criterion

prompt = 'Do you want the animation to be saved? Y/N: ';
str = input(prompt,'s');
if (str == 'Y') || (str == 'y')
    video_par = VideoWriter('outputs\Animation_DLCA.avi'); 
    % Initializing the video file
    video_par.FrameRate = 5;  % Setting frame rate
    open(video_par);  % Opening video file
end

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]);  % Initializing textbar
UTILS.TEXTBAR([1, k_max]);  % Indicating start of marching

for k = 2 : k_max
    
    [par, delt] = MOV.MARCH(par, fl);  % Solving equation of motion
    par.r = MOV.PBC(dom_size, par.r);  % Applying periodic...
    % ...boundary conditions
%     par = COL.GROW(par);  % Checking for collisions and updating particle
%     % structures upon new clusterations
    
    if mod(k-2,t_plt) == 0
        UTILS.PLOTPAR(dom_size, par, 0);  % Plotting every t_plt time steps
        drawnow;  % Drawing the plot at the desired time steps
        pause(0.1);  % Slowing down the animation speed
        if (str == 'Y') || (str == 'y')
            frame_now = getframe(fig_anim);  %  Capturing current frame
            writeVideo(video_par, frame_now);  %  Saving the video
        end
    end
    
    time(k) = time(k-1) + delt; % Updating time
    
    UTILS.TEXTBAR([k, k_max]);  % Updating textbar
    
end

if (str == 'Y') || (str == 'y')
    close(video_par);  % Closing the video file
end


