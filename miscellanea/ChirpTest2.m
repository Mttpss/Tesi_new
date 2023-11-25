clear all
fc = 140e6;
phi0 = 0;
theta0 = 0;
T = 1e-5;
B = 300e6;
f0 = fc;
f1 = fc + B;

M = 5;
Tsw = 2e-5;
fs = fc*10;
Ts = 1/fs;
NT = ceil(T/Ts);
NTsw = ceil(Tsw/Ts);
NTdown = ceil((Tsw-T)/(4*Ts));

t = 0:Ts:(M*NTsw*Ts-Ts);

cUp = (f1-f0)/(NT*Ts);
cDown = (f0-f1)/(NTdown*Ts);

t_up = 0:Ts:((NT-1)*Ts);
t_down = 0:Ts:((NTdown-1)*Ts);

% x = zeros(NTsw*M,1);
% for l = 1:M
%     x(((l-1)*NTsw+1):(NT+(l-1)*NTsw)) = sin(2*pi*(f0*t(((l-1)*NTsw+1):(NT+(l-1)*NTsw))+cUp/2*(t_up.^2)));
%     x((NT+1+(l-1)*NTsw):(NT+NTdown+(l-1)*NTsw)) = sin(2*pi*(f1*t((NT+1+(l-1)*NTsw):(NT+NTdown+(l-1)*NTsw))+cDown/2*(t_down.^2)));
%     x((NT+NTdown+1+(l-1)*NTsw):(NTsw*l)) = zeros(NTsw-NT-NTdown,1);
% end

x = zeros(NTsw*M,1);
for l = 1:M
    x(((l-1)*NTsw+1):(NT+(l-1)*NTsw)) = exp(1j*2*pi*(f0*t(((l-1)*NTsw+1):(NT+(l-1)*NTsw))+cUp/2*(t_up.^2)));
    x((NT+1+(l-1)*NTsw):(NTsw*l)) = zeros(NTsw-NT,1);
end

% y = zeros(NTsw*M,1);
% for l = 1:M
%     y(((l-1)*NTsw+1):(NT+(l-1)*NTsw)) = exp(1j*2*pi*(f0*(t(((l-1)*NTsw+1):(NT+(l-1)*NTsw))-delay(((l-1)*NTsw+1):(NT+(l-1)*NTsw)))+cUp/2*(t_up.^2)));
%     y((NT+1+(l-1)*NTsw):(NTsw*l)) = zeros(NTsw-NT,1);
% end

% x_linear = reshape(x,NTsw*M,1);
x_phase = angle(x);
x_freq = diff(x_phase)/(Ts*2*pi);
hold on
% spectrogram(x,100,80,100,1/Ts,'yaxis')
% spectrogram(y,100,80,100,1/Ts,'yaxis')

% plot(t(1:(end-1)),x_freq)