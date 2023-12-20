function [y,x] = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M)
%ECHOSINGLETARGET Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;
phi0 = 0;
theta0 = 0;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);
t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = 0:Ts:((NT-1)*Ts);
y = zeros(size(t));
d_t = [(d(1) + v(1)*t)' (d(2) + v(2)*t)'];
L = sqrt(10^(-fspl(2*d,lambda_c)/20)); % Free space path loss of the echo
R = dis + v*t;


x = zeros(size(t));
for l = 1:M
    x = 0;
end

y = zeros(NT,M); % Beat signal simulation
for l = 1:M
    y(:,l) = L.*exp(1j*2*pi*(2*fc*d/c + 2*fc*l*NTsw*Ts*v/c + (2*fc*v/c + 2*m*d/c + 2*m*l*NTsw*Ts*v/c)*t_up + (2*m*v/c)*(t_up.^2)));
end

end

