% "main_test" is a script to examine and monitor the performance and...
%   ...outputs of "MCEM" program.

clc
clear
clf('reset')
close all

%% Part 1. Setting the initial domain conditions

[params_ud, params_const] = TRANSP.INIT_PARAMS(); % Initializing the...
    % ...physical parameters to be used in the simulations

[pars, fl] = TRANSP.INIT_DOM(params_ud, params_const); % Initializing...
    % ...the fluid and particle transport variables

% Calculating the primary particle size and number distributions
[pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(4), params_ud.Value(5:6),...
    params_ud.Value(7:9));

% Initializing the primary particle field; Assigning the indices and sizes
pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3)], pars.n);

% Generating the initial aggregates (if applicable)
pars.pp = PAR.INIT_MORPH(pars.pp);

% Assigning the particle initial locations
pars = PAR.INIT_LOC(params_ud.Value(1:3), pars);

% Finding the equivalent particle sizes
pars = PAR.SIZING(pars);

% Finding the initial mobility propeties of particles
pars = TRANSP.MOBIL(pars, fl, params_const);

% Assigning the particle initial velocities
pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const);

% Finding the initial nearest neighbors
ind_trg = (1 : params_ud.Value(4))'; % Indicices of target particles
coef_trg = 5 .* ones(params_ud.Value(4), 1); % Neighboring enlargement...
    % ...coefficients
pars.nnl = COL.NNS(pars, ind_trg, coef_trg);

% % Converting the global structure to aggregate objects
% pars = PAR.PAR2AGG(pars);

disp("The computational domain is successfully initialized...")

% Visualizing the initial particle locations and velocities, and nearest...
    % ...neighbor lists
% figure
% h0_pose = UTILS.PLOTPARS(pars, params_ud.Value(1:3),...
%     'equivalent_volumetric_size', 'on', 'velocity_vector', 'on',...
%     'render', 'on');
% 
% figure
% ind_trg_test = (randperm(params_ud.Value(4),...
%     min([params_ud.Value(4), 2])))'; % A random portion of target...
%         ...indices for nearest neighbor testing
% h0_nntest = UTILS.PLOTNN(pars, params_ud.Value(1:3), ind_trg_test,...
%     coef_trg(ind_trg_test));
% 
% figure
% h0_3d = UTILS.RENDER(pars);

%% Part 2: Simulating the particle aggregations

k_max = 200; % Marching index limit
time = zeros(k_max,1);
% t_plt = 10; % Defining a plotting timeframe criterion
% t_nns = 10; % The timeframe for nearest neighbor search

% prompt = 'Do you want the aggregation animation to be saved? Y/N: ';
% str = input(prompt,'s'); % Animation saving variable (yes/no)
% str = 'y';
% if (str == 'Y') || (str == 'y')
%    % Checking the existence of animation directory
%     if exist('outputs', 'dir') ~= 7
%         mkdir('outputs');
%     end
%     mkdir('outputs');
%     video_par = VideoWriter('outputs\DLCA_anim.avi'); 
%         % Initializing the video file
%     video_par.FrameRate = 3; % Setting frame rate
%     video_par.Quality = 100; % Setting quality
%     open(video_par); % Opening video file
% elseif (str ~= 'N') && (str ~= 'n')
%     error('Error saving the animation (invalid user response)!\n')
% end

% % Initializing the animation plot
% figure
% h_anim = UTILS.PLOTPARS(pars, params_ud.Value(1:3));

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]); % Initializing textbar
UTILS.TEXTBAR([1, k_max]); % Indicating start of marching

for k = 2 : k_max
    if length(cat(1, pars.n)) > (params_ud.Value(4) / 3) % Checking if...
            % ...the number of aggregates within the domain is reasonable
        
        [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solving...
            % ...equation of motions

        pars = TRANSP.PBC(params_ud.Value(1:3), pars); % Applying...
            % ...periodic boundary conditions

        pars = COL.GROW(pars); % Checking for collisions and updating...
          % ...particle structures upon new clusterations

        pars = PAR.SIZING(pars);
            % Updating the size-related properties

        pars = TRANSP.MOBIL(pars, fl, params_const); % Updating the...
            % ...mobility propeties

    %     if mod(k-2, t_nns) == 0
    %         par.nnl = COL.NNS(par, ind_trg, coef_trg); % Updating the...
    %             % ...nearest neighbors
    %     end

    %     if mod(k-2, t_plt) == 0
    %         h_anim = UTILS.PLOTPARS(pars, params_ud.Value(1:3)); % Plotting...
    %             % ...every t_plt time steps
    %         drawnow; % Drawing the plot at the desired time steps
    %         pause(0.1); % Slowing down the animation speed
    %         if (str == 'Y') || (str == 'y')
    %             framenow = getframe(h_anim, [0, 0, 1000, 892.1]);
    %                 % Capturing current frame
    %             writeVideo(video_par, framenow); % Saving the video
    %         end
    %     end

        time(k) = time(k-1) + delt; % Updating time

        UTILS.TEXTBAR([k, k_max]); % Updating textbar
    
    end
end

% if (str == 'Y') || (str == 'y')
%     close(video_par); % Closing the video file
% end

%% Part 3: Postprocessing the results

% Obtaining fractal properties
[df, kf] = PAR.FRACTALITY(pars);

figure
h_3d = UTILS.RENDER(pars); % Morphology of the final population
