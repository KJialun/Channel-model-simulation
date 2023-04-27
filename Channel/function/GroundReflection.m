function [E,reflection_distance,transmit_gain,receive_gain] = GroundReflection(distance)
frequency_carrier = 26*10^9;   % 26 GHz
light_speed = 3*10^8;
htx=2;
hrx=2;
reflection_distance=sqrt((htx+hrx)^2+distance^2);
Theta=(pi/2)-atan((htx+hrx)/distance);
Theta_Gnd=pi-Theta;
transmit_gain=(16/(3*pi))*(sin(Theta_Gnd))^3;
receive_gain=(16/(3*pi))*(sin(Theta_Gnd))^3;
refletion_factor= (cos(Theta)-0.5*sqrt(1-0.25*(sin(Theta)*sin(Theta))))/(cos(Theta)+0.5*sqrt(1-0.25*(sin(Theta)*sin(Theta))));
E=sqrt(transmit_gain*receive_gain)*(refletion_factor*(exp((-1j*2*pi*frequency_carrier/light_speed)*reflection_distance)/reflection_distance));
end

