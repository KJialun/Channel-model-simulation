function [Theta,refletion_factor] = reflection_factor(A,B,C)
%REFLECTION_FACTOR 此处显示有关此函数的摘要
%   此处显示详细说明
v_1 = [A,0] - [B,0];
v_2 = [C,0] - [B,0];
Theta = rad2deg(atan2(norm(cross(v_1, v_2)), dot(v_1, v_2)));
if Theta > 90
    Theta=180-Theta;
end

Theta =  90 - Theta;

refletion_factor= (cos(deg2rad(Theta))-2*sqrt(1-0.25*(sin(deg2rad(Theta))*sin(deg2rad(Theta)))))/(cos(deg2rad(Theta))+2*sqrt(1-0.25*(sin(deg2rad(Theta))*sin(deg2rad(Theta)))));


end

