clear all
load('FMCWParameterSetup_1.mat');
c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);
t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = 0:Ts:((NT-1)*Ts);
x = zeros(NTsw*M,1);
for l = 1:M
    x(((l-1)*NTsw+1):(NT+(l-1)*NTsw)) = exp(1j*2*pi*(fc*t(((l-1)*NTsw+1):(NT+(l-1)*NTsw))+m/2*(t_up.^2)));
    x((NT+1+(l-1)*NTsw):(NTsw*l)) = zeros(NTsw-NT,1);
end

Bi = (B-50*B/100):B/100:(B+50*B/100);
Ri = zeros(size(Bi));
peaki = zeros(size(Bi));
x_i = zeros(NT,M);
y_i = zeros(NT,M);
d = 200;
v = -15;
y_tf_ts = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M);
y_d_v = zeros(N,M);
y_d_v_sample = zeros(N,M,5);
Bi_sample = zeros(1,5);
i = 1;
x_tot = x;
for n = 1:size(Bi,2)
    disp('Siamo alla iterazione: ');
    disp(n);
    fci = fc;
    di = 100;
    vi = 30;
    [y_i,x_i] = echoSingleInterference(di,vi,fc,fci,fs,Bi(n),NT,NTsw,M);
   
    Ri(n) = max(abs(xcorr(x,x_i,'normalized')),[],'all');
    y_tf_ts = y_tf_ts + y_i;
    y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);
    peaki(n) = max(abs(y_d_v),[],'all');
    if (i-1)*(101-1)/4+1 == n
        y_d_v_sample(:,:,i) = y_d_v;
        Bi_sample(i) = Bi(n);
        x_tot = x_tot + x_i;
        i = i+1;
    end
    % desiredSNR = 0;
    % noise = awgn(ones(1,NT*M),0) - 1;
    % noise = lowpass(noise,B/fs);
    % noise = reshape(noise,NT,M);  
end
spectrogram(x_tot,1024*8,6000,1024*8,1/Ts,'yaxis')
v_ax = ((-M/2):(M/2))*delta_v;
d_ax = (0:(N-1))*delta_d;
for i = 1:5
    figure
    imagesc(v_ax(round(M/4):round(3*M/4)),...
    flipud(d_ax(1:round(size(d_ax,2)/2))'),...
    flipud(abs(y_d_v_sample((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4),i)/...
    max(abs(y_d_v_sample((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4),i)),[],'all'))))
    xlabel('Velocity [m/s]')
    ylabel('Distance [m]')
    cb = colorbar;
    % cb.Label.String = ('[dB]');
    set(gca,'YDir','normal')
    title(['Banda = ', num2str(Bi_sample(i))], 'FontSize', 14);
end

figure
tiledlayout(2,1)

% Top plot
ax1 = nexttile;
plot(ax1,Bi,Ri)
title(ax1,'Banda - Correlazione')
ax1.FontSize = 14;
ax1.XColor = 'red';

% Bottom plot
ax2 = nexttile;
plot(ax2,Bi,peaki)
title(ax2,'Banda - Picco')
ax2.FontSize = 14;
ax2.XColor = 'blue';




