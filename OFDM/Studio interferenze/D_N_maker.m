function [D_N] = D_N_maker(x,N)
%D_N_MAKER Summary of this function goes here
%   Detailed explanation goes here
d_N = zeros(1,N);
for i = 1:N
    d_N(i) = exp(1i*2*pi*(x)*(i-1));
end
D_N = diag(d_N);
end

