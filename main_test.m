% "main_test" is a script to examine and monitor the performance and...
    % ...outputs of "MCEM" program.

clc
clear
clf('reset')
close all

%% Part 1: Defining the constant physical parameters of the problem

% Creating the table of constant parameters
Name = {'rho_bc'; 'M_air'; 'kb'; 'Na'; 'Ru'};
Value = [1.8e3; 28.97e-3; 1.381e-23; 6.022e23; 8.314];
Unit = {'kg/m3'; 'kg/mol'; 'j/k'; 'mol^-1'; 'j/mol.k'};
Description = {'Black Carbon bulk density'; 'Air molar mass';...
    'Boltzmann constant'; 'Avogadro constant'; 'Universal gas constant'};
params_const = table(Name, Value, Unit, Description);

clear Name Value Unit Description % Cleaning the worksapce

%% Part 2: Importing the user-defined parameters of domain, particle and flow

% Extracting data from the input file
file_id = fopen('inputs\MCEM_Params.txt','r');
import_data = textscan(file_id,'%s%f%s%s','Headerlines',2,...
    'Delimiter','\t','MultipleDelimsAsOne',1);
fclose(file_id);

% Creating the table of user-defined parameters
Name = import_data{1,1};
Value = import_data{1,2};
Unit = import_data{1,3};
Description = import_data{1,4};
params_ud = table(Name, Value, Unit, Description);
% NOTE: See the "Description" column for more info on these parameters.

clear file_id import_data ans Name Value Unit Description % Cleaning the...
    % ...worksapce

%% Part 3: Initializing the fluid and particle transport variables

% Declaring the fluid structure and loading the corresponding imported...
    % ...properties
fl = struct('temp', params_ud.Value(10), 'v', params_ud.Value(11:3),...
    'p', params_ud.Value(14), 'mu', [],'lambda', []);
% "fl" is a structure containing the main physical properties of the...
    % ...background fluid. It is used to consider and study the...
    % ...particle-fluid interactions.
% Fields:
    % temp: Fluid temperature
    % v: ~ velocity
    % p: ~ pressure
    % mu: ~ viscosity
    % lambda: ~ mean free path

[fl.mu, fl.lambda] = TRANSP.FLPROPS(fl,params_const); % Calculating the...
    % ...fluid viscosity and mean free path

% Declaring the particle structure
pars = struct('pp', [], 'n', [], 'd', [], 'r', [], 'v', [], 'm', [],...
    'rho', [], 'delt', [], 'tau', [], 'f', [], 'diff', [],...
    'lambda', [], 'kn', [], 'nnl', []);
% "par" is a structure conating the main physical properties of...
    % ...independent particles (whether being single monomers or...
    % ...aggregates) as well as the properties of their constituent...
    % ...primary particles (if applicable).
% Fields:
    % pp: Primary particle characteristics (an N*5 matrix containing...
        % ...their index, size, and location)
    % n: Number of primaries within the particles
    % d: Independent particles equivalent diameter
    % r: ~ spatial location
    % v: ~ velocity
    % m: ~ mass
    % rho: ~ effective density
    % delt: ~ motion time-step
    % tau: ~ relaxation time
    % f: ~ friction factor
    % diff: ~ diffusivity
    % lambda: ~ diffusive mean free path
    % kn: Knudsen number (both kinetic and diffusive)
    % nnl: ~ nearest neighbor list
% NOTE: The rows of each field correspond to characteristics of each...
    % ...independent particle.

% Calculating the primary particle size and number distributions
[pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(4), params_ud.Value(5:6),...
    params_ud.Value(7:9));

% Initializing the primary particle field; Assigning the indices and sizes
pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3)], pars.n);

% Generating the initial aggregates (if applicable)
[pars.pp, pars.d] = PAR.INIT_MORPH(pars.pp);

% Assigning the particle initial locations
pars = PAR.INIT_LOC(params_ud.Value(1:3), pars);

% % Finding the equivalent particle sizes
% par.d = TRANSP.SIZING(par.pp, par.n);

% Finding the initial mobility propeties of particles
pars = TRANSP.MOBIL(pars, fl, params_const);

% Assigning the particle initial velocities
pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const);

% Finding the initial nearest neighbors
ind_trg = (1 : params_ud.Value(4))'; % Indicices of target particles
coef_trg = 5 .* ones(params_ud.Value(4), 1); % Neighboring enlargement...
    % ...coefficients
pars.nnl = COL.NNS(pars, ind_trg, coef_trg);

% pars = PAR.PAR2AGG(pars);  % Converting to aggregate objects

disp("The computational domain is successfully initialized...")

% Visualizing the initial particle locations and velocities, and nearest...
    % ...neighbor lists
% figure
% fig_parinit = UTILS.PLOTPAR(params_ud.Value(1:3), pars,...
%     'equivalent_volumetric_size', 'on', 'velocity_vector', 'on');

% figure
% ind_trg_test = (randperm(params_ud.Value(4),...
%     min([params_ud.Value(4), 5])))'; % A random portion of target...
%         % ...indices for nearest neighbor testing
% [fig_parinit, fig_nntest] = UTILS.PLOTPAR(params_ud.Value(1:3), pars,...
%     'equivalent_volumetric_size', 'on', 'velocity_vector', 'on',...
%     'nearest_neighbor', 'on', 'target_index', ind_trg_test,...
%     'target_coefficient', coef_trg(ind_trg_test));

% figure
% fig_3d = UTILS.RENDER(params_ud.Value(1:3), pars.pp, pars.n);

%% Part 4: Simulating the particle aggregations

k_max = 50; % Marching index limit
time = zeros(k_max,1);
t_plt = 1; % Defining a plotting timeframe criterion
t_nns = 10; % The timeframe for nearest neighbor search

% prompt = 'Do you want the aggregation animation to be saved? Y/N: ';
% str = input(prompt,'s'); % Animation saving variable (yes/no)
str = 'y';
if (str == 'Y') || (str == 'y')
%    % Checking the existence of animation directory
%     if exist('outputs', 'dir') ~= 7
%         mkdir('outputs');
%     end
    mkdir('outputs');
    video_par = VideoWriter('outputs\DLCA_anim.avi'); 
        % Initializing the video file
    video_par.FrameRate = 3; % Setting frame rate
    video_par.Quality = 100; % Setting quality
    open(video_par); % Opening video file
end

% Initializing the animation
figure
fig_anim = UTILS.PLOTPAR(params_ud.Value(1:3), pars);

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]); % Initializing textbar
UTILS.TEXTBAR([1, k_max]); % Indicating start of marching

for k = 2 : k_max
    
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solving...
        % ...equation of motions
    
    pars = TRANSP.PBC(params_ud.Value(1:3), pars); % Applying periodic...
        % ...boundary conditions
    
    pars = COL.GROW(pars); % Checking for collisions and updating...
      % ...particle structures upon new clusterations
    
%     par.d = TRANSP.SIZING(par.pp, par.n); % Updating the equivalent sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Updating the mobility...
        % ...propeties
    
%     if mod(k-2, t_nns) == 0
%         par.nnl = COL.NNS(par, ind_trg, coef_trg); % Updating the...
%             % ...nearest neighbors
%     end

    if mod(k-2, t_plt) == 0
        fig_anim = UTILS.PLOTPAR(params_ud.Value(1:3), pars); % Plotting...
            % ...every t_plt time steps
        drawnow; % Drawing the plot at the desired time steps
        pause(0.1); % Slowing down the animation speed
        if (str == 'Y') || (str == 'y')
            frame_now = getframe(fig_anim, [0, 0, 1000, 892.1]);
                % Capturing current frame
            writeVideo(video_par, frame_now); % Saving the video
        end
    end
    
    time(k) = time(k-1) + delt; % Updating time
    
    UTILS.TEXTBAR([k, k_max]); % Updating textbar
    
end

if (str == 'Y') || (str == 'y')
    close(video_par); % Closing the video file
end
