function A = delta(T,T_N)
%DELTA Summary of this function goes here
%   Detailed explanation goes here

A=dirac(round(T-T_N,10));
A(A==inf)=1;
end

