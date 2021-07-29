addpath('Functions/')

clear all; 
close all;

num_of_mol = 1;

total_time = 40; % time in minutes
time_step = 1; % in seconds
t = 0:time_step:total_time*60; % in seconds

% Starting time to each cluster is Normal distribution. Values taken from
% fit to real data
mu_start_time = 13; % in minutes
sigma_start_time = 3; % in minutes

% Create starting time dist (and sample from it)
pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / time_step); % notice conversion to index



z = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

figure 
hold on

for i=1:10
    
    zeros_over_time = startingHeightTest(z(i), t, sigmoid_start, time_step);
    plot(t, zeros_over_time, 'LineWidth', 1, 'DisplayName', sprintf('z= %0.1f', z(i)));
end

hold off
legend('show');
title('Number of sensors with protons in top 25% of all');
%savefig('instant_nH_zeros_no_diff.fig');



