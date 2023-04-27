% Set user-defined simulation parameters.
% Speed of light:
light_speed = 3*10^8;
% Distance between points along ray:
ray_resolution = 0.1;
% Maximum number of ray collisions:
ray_collisions = 1;
% Carrier frequency:
frequency_carrier = 26*10^9;   % 26 GHz
EIRP=0.25;
% Transmitter gain:
transmit_gain = (16/(3*pi))*(sin(pi/2))^3;            
P_TX=EIRP/transmit_gain;
% Receiver gain:
receive_gain = (16/(3*pi))*(sin(pi/2))^3;              

red = [1,10/51,10/51];
blue = [14/51,32/51,1];
green = [0,40/51,10/51];
yellow = [50/51,40/51,4/17];