function [y,x] = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M)
%ECHOSINGLETARGET Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);

t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = (0:Ts:((NT-1)*Ts))';
R = d + v*t;
L = (sqrt(10.^(-fspl(2*R,lambda_c)/20)))';
L = reshape(L,NTsw,M);

x = zeros(NTsw,M);
f = zeros(NTsw,M);
for l = 1:M
    x(1:NT,l) = exp(1i*2*pi*(fc*t_up + m/2*((t_up).^2)));
    x((NT+1):end,l) = zeros(NTsw-NT,1);
end
x = reshape(x,NTsw*M,1);

y = zeros(NT,M); % Beat signal simulation
for l = 1:M
    y(:,l) = L(1:NT,l).*exp(1j*2*pi*(2*fc*d/c + 2*fc*l*NTsw*Ts*v/c + (2*fc*v/c + 2*m*d/c + 2*m*l*NTsw*Ts*v/c)*t_up + (2*m*v/c)*(t_up.^2)));
end

end

