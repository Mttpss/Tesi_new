clear all
close all
load('RadarParameters\FMCWParameterSetup_1.mat');
c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);
t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = 0:Ts:((NT-1)*Ts);
numofpoint = 10; % deve essere pari

%% TARGET'S ECHOES GENERATION

d = 100;
v = 15;
[y_target,x,f] = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M);
sin_t = (x - conj(x))/(1i*2);

%% NOISE GENERATION

desiredSNR = 0;
noise_pre = awgn(y_target,0) - y_target;
noise_pre_filtered = noise_pre;
if B/fs < 1
    noise_pre_filtered = lowpass(noise_pre,B/fs);
end
Pt = sum(y_target.*conj(y_target),"all");
a = sqrt((Pt/(sum(noise_pre_filtered.*conj(noise_pre_filtered),"all")))*10^(-desiredSNR/10));
noise = a*noise_pre_filtered;
noise_fft = fft(noise)/fs;
Pn = sum(noise.*conj(noise),"all");
SNRdB = 10*log10(Pt/Pn);

%% INTERFERENCE GENERATION

n = 1;
di = 100;
vi = 20;
fci = fc - B/2;
Bi = B*3;
Ti = NT*Ts;
NTi = floor(Ti/Ts);
Tswi = NTsw*Ts;
NTswi = floor(Tswi/Ts);
Mi = M;
% offset = 0;
t_i = 0:Ts:(Mi*NTswi*Ts-Ts);

R = zeros(numofpoint,1);
SIRdB = zeros(numofpoint,1);
y_i_tot = zeros(NTsw,M);
xi_sample = zeros(NTsw*M,1);
yi_sample = zeros(NT,M);
y_d_v_sample = zeros(N,M);
offset_sample = 0;
SIRdB_sample = 0;
R_sample = 0;
offset = (-M*NTsw + 1):(M*NTsw - 1);

l=1;
mezzo_flag = 0;
% 1:(2*M*NTsw - 1)
for k = 1:numofpoint
    disp('Iterazione ')
    disp(k)
    u = floor(k*(2*M*NTsw+1)/numofpoint - (2*M*NTsw+1)/(numofpoint*2));
    disp(u)
    disp(offset(u))
    [yi,xi] = echoInterferenceFMCW(n,di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset(u));
    offset_ax(k) = offset(u);
    y_i_noise = yi + noise;
    y_tf_ts = y_target + yi + noise;
    sin_i = (xi - conj(xi))/(1i*2);
    R(k) = max(xcorr(sin_t,sin_i,'normalized'),[],'all');
    Pi = sum(yi.*conj(yi),"all");
    SIRdB(k) = 10*log10(Pt/(Pi + Pn));
    if l*10 == k || k == 1 || u == M*NTsw
        if u == M*NTsw
            mezzo_flag = 1;
        end
        xi_sample = cat(3,xi_sample,xi);
        yi_sample = cat(3,yi_sample,yi);
        y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);
        y_d_v_sample = cat(3,y_d_v_sample,y_d_v);
        offset_sample = [offset_sample offset_ax(k)];
        SIRdB_sample = [SIRdB_sample SIRdB(k)];
        R_sample = [R_sample R(k)]
        l = l+1;
    end
end
xi_sample(:,:,1) = [];
yi_sample(:,:,1) = [];
y_d_v_sample(:,:,1) = [];
offset_sample(1) = [];
SIRdB_sample(1) = [];
R_sample(1) = [];
if mezzo_flag == 0
    [yi,xi] = echoInterferenceFMCW(n,di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,0);
    offset_ax = [offset_ax(1:floor(end/2)) , 0 , offset_ax(ceil(end/2):end)];
    y_i_noise = yi + noise;
    y_tf_ts = y_target + yi + noise;
    sin_i = (xi - conj(xi))/(1i*2);
    R = [R(1:floor(end/2)) ; max(xcorr(sin_t,sin_i,'normalized'),[],'all') ; R(ceil(end/2):end)];
    Pi = sum(yi.*conj(yi),"all");
    SIRdB = [SIRdB(1:floor(end/2)) ; 10*log10(Pt/(Pi + Pn)) ; SIRdB(ceil(end/2):end)];
    xi_sample = cat(3,xi_sample(:,:,1:floor(end/2)) , xi , xi_sample(:,:,ceil(end/2):end));
    yi_sample = cat(3,yi_sample(:,:,1:floor(end/2)) , yi , yi_sample(:,:,ceil(end/2):end));
    y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);
    y_d_v_sample = cat(3,y_d_v_sample(:,:,1:floor(end/2)) , y_d_v , y_d_v_sample(:,:,ceil(end/2):end));
    offset_sample = [offset_sample(1:floor(end/2)) 0 offset_sample(ceil(end/2):end)];
    SIRdB_sample = [SIRdB_sample(1:floor(end/2)) 10*log10(Pt/(Pi + Pn)) SIRdB_sample(ceil(end/2):end)];
    R_sample = [R_sample(1:floor(end/2)) max(xcorr(sin_t,sin_i,'normalized'),[],'all')  R_sample(ceil(end/2):end)];
    l = l+1;
end  

%% PLOTTING VARIOUS RESULTS
figure
tiledlayout(2,1)

% Top plot
ax1 = nexttile;
plot(ax1,offset_ax,R)
title(ax1,'Offset - Correlazione')
ax1.FontSize = 14;
ax1.XColor = 'red';

% Bottom plot
ax2 = nexttile;
plot(ax2,offset_ax,SIRdB)
title(ax2,'Offset - SIR')
ax2.FontSize = 14;
ax2.XColor = 'blue';

for k = 1:size(xi_sample,3)
    v_ax = ((-M/2):(M/2))*delta_v;
    d_ax = (0:(N-1))*delta_d;

    figure
    imagesc(v_ax(round(M/4):round(3*M/4)),...
    flipud(d_ax(1:round(size(d_ax,2)/2))'),...
    flipud(20*log10(abs(y_d_v_sample((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4),k)/...
    max(abs(y_d_v_sample((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4),k)),[],'all')))))
    xlabel('Velocity [m/s]')
    ylabel('Distance [m]')
    cb = colorbar;
    cb.Label.String = ('[dB]');
    set(gca,'YDir','normal')
title(['Offset = ', num2str(offset_sample(k)) ,' | SIR = ', num2str(SIRdB_sample(k)), ' dB | Correlation = ', num2str(R_sample(k)) ], 'FontSize', 14);    
    figure
    imagesc(v_ax(round(M/4):round(3*M/4)),...
    flipud(d_ax(1:round(size(d_ax,2)/2))'),...
    flipud(abs(y_d_v_sample((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4),k)/...
    max(abs(y_d_v_sample((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4),k)),[],'all'))))
    xlabel('Velocity [m/s]')
    ylabel('Distance [m]')
    cb = colorbar;
    set(gca,'YDir','normal')
    title(['Offset = ', num2str(offset_sample(k)) ,' | SIR = ', num2str(SIRdB_sample(k)), ' dB | Correlation = ', num2str(R_sample(k)) ], 'FontSize', 14);
end


figure
tiledlayout(1,1)

% Bottom plot
ax3= nexttile;
plot(ax3,R(ceil(end/2):end),SIRdB(ceil(end/2):end))
title(ax3,'Correlation - SIR')
ax3.FontSize = 14;
ax3.XColor = 'blue';




