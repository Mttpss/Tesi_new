function [y,x,f] = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M)
%ECHOSINGLETARGET Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);

t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = (0:Ts:((NT-1)*Ts))';
% y = zeros(size(t));
% L = sqrt(10^(-fspl(2*d,lambda_c)/20)); % Free space path loss of the echo
R = d + v*t;
L = (sqrt(10.^(-fspl(2*R,lambda_c)/20)))';
L = reshape(L,NTsw,M);

% x = zeros(NTsw*M,1);
% for l = 1:M
%     x(((l-1)*NTsw+1):(NT+(l-1)*NTsw)) = exp(1j*2*pi*(fc*t(((l-1)*NTsw+1):(NT+(l-1)*NTsw))+m/2*(t_up.^2)));
%     x((NT+1+(l-1)*NTsw):(NTsw*l)) = zeros(NTsw-NT,1);
% end

x = zeros(NTsw,M);
f = zeros(NTsw,M);
for l = 1:M
    % x(1:NT,l) = exp(1i*2*pi*(fc*(t_up + (l-1)*NTsw*Ts) + m/2*((t_up).^2)));
    x(1:NT,l) = exp(1i*2*pi*(fc*t_up + m/2*((t_up).^2)));
    x((NT+1):end,l) = zeros(NTsw-NT,1);
    f(1:NT,l) = fc + m/2*t_up;
    f((NT+1):end,l) = zeros(NTsw-NT,1);
end
x = reshape(x,NTsw*M,1);
f = reshape(f,NTsw*M,1);


y = zeros(NT,M); % Beat signal simulation
for l = 1:M
    y(:,l) = L(1:NT,l).*exp(1j*2*pi*(2*fc*d/c + 2*fc*l*NTsw*Ts*v/c + (2*fc*v/c + 2*m*d/c + 2*m*l*NTsw*Ts*v/c)*t_up + (2*m*v/c)*(t_up.^2)));
end

end

