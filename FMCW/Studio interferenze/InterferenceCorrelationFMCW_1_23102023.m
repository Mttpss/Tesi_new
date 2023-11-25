clear all
close all
load('FMCWParameterSetup_1.mat');
c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);
t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = 0:Ts:((NT-1)*Ts);

d = 100;
v = 15;
[y_target,x,f] = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M);
x = x.*conj(exp(1i*2*pi*fc*t)');
f = f-fc;
% y_tf_ts = zeros(NT,M);
% spectrogram(x,100,80,1024,1/Ts,'yaxis')
% for l = 1:8
%     figure
%     spectrogram(x(NTsw*(l-1)+1:NTsw*(l-1)+1+NT),100,80,1024,1/Ts,'yaxis')
% end

di = 100;
vi = 20;
fci = fc;
Bi = B;
Ti = NT*Ts;
NTi = floor(Ti/Ts);
Tswi = NTsw*Ts;
NTswi = floor(Tswi/Ts);
Mi = M;
offset = NTsw*M/2;
t_i = 0:Ts:(Mi*NTswi*Ts-Ts);

R = zeros(NTsw*M,1);


for offset = (-M*NTsw):(M*NTsw)
    
end

[yi,xi,fi] = echoSingleInterferenceFMCW(di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset);
fi = fi(1:NTsw*M,1)-fci;;
% xi = xi(1:NT*M,1);
% x_tot = x + x_i;
% spectrogram(x_tot,100,80,1024,1/Ts,'yaxis')

% y_i = zeros(NT,M);

% Ri = max(abs(xcorr(x,x_i,'normalized')),[],'all');
y_tf_ts = y_tf_ts + yi;
y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);
% peak = max(abs(y_d_v),[],'all');

%% ADDING NOISE

desiredSNR = 0;
noise_pre = awgn(y_target,0) - y_target;
noise_pre_filtered = noise_pre;
noise_pre_filtered = lowpass(noise_pre,B/fs);
a = sqrt((sum(y_target.*conj(y_target),2)/(sum(noise_pre_filtered.*conj(noise_pre_filtered),2)))*10^(-desiredSNR/10));
noise = a*noise_pre_filtered;
noise_fft = fft(noise)/fs;

%% PLOTTING VARIOUS RESULTS

v_ax = ((-M/2):(M/2))*delta_v;
d_ax = (0:(N-1))*delta_d;

figure
tiledlayout(2,1)

% Top plot
ax1 = nexttile;
plot(ax1,t,f)
title(ax1,'Frequenza segnale radar')
ax1.FontSize = 14;
ax1.XColor = 'red';

% Bottom plot
ax2 = nexttile;
plot(ax2,t,fi)
title(ax2,'Frequenza interferenza')
ax2.FontSize = 14;
ax2.XColor = 'blue';

figure
imagesc(v_ax(round(M/4):round(3*M/4)),...
flipud(d_ax(1:round(size(d_ax,2)/2))'),...
flipud(20*log10(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))/...
max(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))),[],'all')))))
xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')
title(['Banda = ', num2str(Bi)], 'FontSize', 14);

figure
imagesc(v_ax(round(M/4):round(3*M/4)),...
flipud(d_ax(1:round(size(d_ax,2)/2))'),...
flipud(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))/...
max(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))),[],'all'))))
xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')
title(['Banda = ', num2str(Bi)], 'FontSize', 14);
