function [P] = P_maker(alpha,gamma,Nc,Nofdm)
%P_MAKER Summary of this function goes here
%   Detailed explanation goes here
P = zeros(Nc,Nofdm);
for n = 0:(Nc-1)
    for u = 0:(Nofdm-1)
        P(n+1,u+1) = exp(-1i*2*pi*n*gamma*u*alpha);
    end
end
end

