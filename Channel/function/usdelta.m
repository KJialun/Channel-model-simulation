function A = usdelta(T,T_N)
%DELTA Summary of this function goes here
%   Detailed explanation goes here
A=dirac(round(T-T_N,8));
A(A==inf)=1;
end
