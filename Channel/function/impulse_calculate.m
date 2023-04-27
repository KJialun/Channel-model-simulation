function h_T = impulse_calculate(T,E,T_n)
parameters
received_power=(light_speed/(frequency_carrier*4*pi))^2*P_TX.*abs(E).^2;
a=sqrt(received_power./P_TX);
% a=ones(length(received_power),1);
h_T=sum(a.*exp((-1j*2*pi*frequency_carrier).*T_n).*delta(T,T_n));
h_T=abs(h_T);
end

