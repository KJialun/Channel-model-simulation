function h_TDL = US_TDL_calc(T,E,T_n,bandwidth)
parameters

delta_T=1/(2*bandwidth);
max_delay=max(T_n);
L=round(max_delay/delta_T);

received_power=(light_speed/(frequency_carrier*4*pi))^2*P_TX.*abs(E).^2;
a=sqrt(received_power./P_TX);

h_TDL=0;

for i=0:L
h_L=sum(a.*exp((-1j*2*pi*frequency_carrier).*T_n).*usdelta(T_n,i*delta_T));
h_TDL=h_TDL+h_L*delta(T,i*delta_T);
end

h_TDL=abs(h_TDL);


end

