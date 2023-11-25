function [F_N] = F_N_matrix_maker(N)
%FOURIER_MATRIX_MAKER Summary of this function goes here
%   Detailed explanation goes here
F_N = zeros(N);
for k = 0:(N-1)
    for n = 0:(N-1)
        F_N(k+1,n+1) = (1/sqrt(N))*exp(-1i*2*pi*k*n/N);
    end
end
end

