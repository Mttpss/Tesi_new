function [y,x] = echoSingleTarget(d,v,Pt,fc,fs,B,NT,NTsw,M)
%ECHOSINGLETARGET Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);

t = (0:Ts:(M*NTsw*Ts-Ts))';
t_up = (0:Ts:((NT-1)*Ts))';
d_2d = [d(1) + v(1)*t , d(2) + v(2)*t];
D = sqrt(d_2d(:,1).^2 + d_2d(:,2).^2);
V = v(1)*(d_2d(:,1)./D) + v(2)*(d_2d(:,2)./D);
L = (sqrt(10.^(-fspl(2*D,lambda_c)/20)))';
D = reshape(D,NTsw,M);
L = reshape(L,NTsw,M);

x = zeros(NTsw,M);
for l = 1:M
    x(1:NT,l) = exp(1i*2*pi*(fc*t_up + m/2*((t_up).^2)));
    x((NT+1):end,l) = zeros(NTsw-NT,1);
end
x = reshape(x,NTsw*M,1);
a = sqrt((Pt/(sum(x.*conj(x),"all"))));

y = zeros(NT,M); % Beat signal simulation
for l = 1:M
    y(:,l) = L(1:NT,l).*exp(1i*2*pi*(2*fc*D(1:NT,l)/c + 2*fc*l*NTsw*Ts*V(1)/c + (2*fc*V(1)/c + 2*m*D(1:NT,l)/c + 2*m*l*NTsw*Ts*V(1)/c)*t_up + (2*m*V(1)/c)*(t_up.^2)));
end
y = a*y;
end

