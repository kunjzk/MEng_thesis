%% Clear workspace and load functions
clear all; 
close all;
addpath('./Functions/');

%% Initialise constants

PLOT = 1; % flag to plot figures or not

N = 32; % N x N array
num_of_mol = 20; % # of initial molecules
sensor_dim = 4E-3; % in meters
solution_height = 2.5E-3; % in meters
M = 32; % resolution of Z axis. I.e. N x N x M space
speed_of_front = 1.73E-6 ; % meters per second

total_time = 40; % time in minutes
time_step = 5; % in seconds
t = 0:time_step:total_time*60; % in seconds

% Starting time to each cluster is Normal distribution...
mu_start_time = 13; % in minutes
sigma_start_time = 3; % in minutes

% Number of protons released is proportional to area of propagation front.
% 'burst' is the proportionality constant..
burst = 50;

est_tot_protons = estimate_total_protons(N, M, sensor_dim, solution_height, speed_of_front, t, burst);

% Create starting time dist (and sample from it)
pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / time_step); % notice conversion to index 

% Initialise those "clusters"
if exist('clusters_init','var')
    clear clusters_init
end
clusters_init(1, num_of_mol) = struct();
for i = 1:num_of_mol
    clusters_init(i).centre_x = 1 + N*rand;
    clusters_init(i).centre_y = 1 + N*rand;
    clusters_init(i).centre_z = 1 + M*rand;
    clusters_init(i).radius_h = 0.01;
    clusters_init(i).radius_v = BuMeters(BuMeters(0.01,0,sensor_dim,N), 1,solution_height,M);
    clusters_init(i).protons_x = [];
    clusters_init(i).protons_y = [];
    clusters_init(i).protons_z = [];
end

%% Visualise

if(PLOT==1)
    figure;
    hh = scatter3(cat(1,clusters_init.centre_x),...
                  cat(1,clusters_init.centre_y),...
                  cat(1,clusters_init.centre_z), 'LineWidth', 3);
              
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([1 N]); ylim([1 N]); zlim([1 M]);
    title({'Starting Position of Initial Molecules in 3D Space';...
          ['Number of Molecules = ' num2str(num_of_mol)]}, 'FontSize', 20);
    
    % Aesthetics...
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;
    
    % Plot something that looks like an ISFET array :P
    hold on; 
    [X,Y] = meshgrid(1:N);
    h = surf(X, Y, ones(N));
    set(h,'facealpha',0.5);
end