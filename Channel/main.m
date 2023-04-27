% Close all open figure windows:
addpath(genpath('function'));
close all
% Remove existing variables from memory:
clear
% Refresh command window:
% clc
transmit_pos = [35 ,95];
receive_pos = [70,10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters
environment
draw_rays_on_map
channel_impulse_response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_1d
% whole_map_calc  %% calculate P_RX at all of receivers 
% plot_heatmap