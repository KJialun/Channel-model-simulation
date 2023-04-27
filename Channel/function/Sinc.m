function OUT = Sinc(T_N,delta_T,L,f_sampling)
%SINC Summary of this function goes here
%   Detailed explanation goes here
delta_T=round(delta_T.*f_sampling);
OUT=sinc(2*(T_N-L*delta_T));
end

