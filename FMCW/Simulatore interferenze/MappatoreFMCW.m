clear all
close all
%% SETTING FLAG FOR DESIRED MAPS
linear_map_flag = 0;
dB_map_flag = 1;

%% SETTING RADAR PARAMETERS

loadRadarSetup_flag = 0;
if loadRadarSetup_flag ~= 0
    load('RadarParameters\FMCWParameterSetup_1.mat');
else
    du_ideal = 300;
    vu_ideal = 95;
    delta_d_ideal = 0.3/3.6;
    delta_v_ideal = 1/3.6;
    B = 500e6;
    fc = 140e9;
    fs = 1000e6;

    [phyReal_flag,du,vu,delta_d,delta_v,Ts,Tsw,NTsw,T,NT,N,M] = ...
        parametrizationFMCW(du_ideal,vu_ideal,delta_d_ideal,delta_v_ideal,B,fc,fs);

    if phyReal_flag == 0
        disp(['Non è possibile raggiungere contemporaneamente i requisiti' ...
            ' di massima distanza e velocità non-ambigue']);
        return;
    end

end

c = 3e8;
lambda_c = c/fc;
Ts = 1/fs;
m = B/(NT*Ts);
t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = 0:Ts:((NT-1)*Ts);

%% TARGET'S ECHOES GENERATION

d = 100;
v = 15;
[y_target,x] = echoSingleTarget(d,v,fc,fs,B,NT,NTsw,M);
sin_t = (x - conj(x))/(1i*2);

%% NOISE GENERATION

desiredSNR = 0;
noise_pre = awgn(y_target,0) - y_target;
noise_pre_filtered = noise_pre;
% if B/fs < 1
%     noise_pre_filtered = lowpass(noise_pre,B/fs);
% end
Pt = sum(y_target.*conj(y_target),"all");
a = sqrt((Pt/(sum(noise_pre_filtered.*conj(noise_pre_filtered),"all")))*10^(-desiredSNR/10));
noise = a*noise_pre_filtered;
noise_fft = fft(noise)/fs;
Pn = sum(noise.*conj(noise),"all");
SNRdB = 10*log10(Pt/Pn);

%% INTERFERENCE POINT SETUP

load_point_flag = 0;

if load_point_flag == 1
    load('point.mat');
else
    numofpoint = 20; % must be even
    fci_flag = 0;
    Bi_flag = 0;
    Mi_flag = 0;
    Ti_flag = 0;
    Tswi_flag = 0;
    offset_flag = 1;

    if Bi_flag ~= 0
        Bi = linspace(0,B,numofpoint/2);
        Bi(end) = [];
        Bi = [Bi linspace(B,B*2,numofpoint/2)];
    else
        Bi = ones(1,numofpoint-1)*B;
    end

    if fci_flag ~= 0
        fci = linspace(fc-B,fc,numofpoint/2);
        fci(end) = [];
        fci = [fci linspace(fc,fc+B,numofpoint/2)];
    else
        fci = ones(1,numofpoint-1)*fc;
    end

    if Mi_flag ~= 0
        Mi = floor(linspace(0,M,numofpoint/2));
        Mi(end) = [];
        Mi = [Mi floor(linspace(M,2*M,numofpoint/2))];
    else
        Mi = ones(1,numofpoint-1)*M;
    end

    if Ti_flag ~= 0
        Ti = linspace(NT*Ts/10,NT*Ts,numofpoint/2);
        Ti(end) = [];
        Ti = [Ti linspace(NT*Ts,NT*Ts*2,numofpoint/2)];
    else
        Ti = ones(1,numofpoint-1)*NT*Ts;
    end

    if Tswi_flag ~= 0
        Tswi = floor(linspace(0,(NTsw-NT),numofpoint/2))*Ts;
        Tswi(end) = [];
        Tswi = [Tswi floor(linspace((NTsw-NT),(NTsw-NT)*2,numofpoint/2))*Ts];
        Tswi = Ti + Tswi;
    else
        Tswi = ones(1,numofpoint-1)*(NTsw-NT)*Ts+Ti;
    end

    if offset_flag ~= 0
        offset = floor(linspace(-floor(max(Tswi,[],"all")/Ts*max(Mi,[],"all")),0,numofpoint/2));
        offset(end) = [];
        offset = [offset floor(linspace(0,NTsw*M,numofpoint/2))];
    else
        offset = zeros(1,numofpoint-1);
    end
end

%% INTERFERENCE GENERATION

n = 1;
di = 100;
vi = 20;

R = zeros(numofpoint-1,1);
SIR = zeros(numofpoint-1,1);
SIRdB = zeros(numofpoint-1,1);

for k = 1:(numofpoint-1)
    disp('Iterazione ')
    disp(k)
    [yi,xi] = echoInterferenceFMCW(x,Pt,n,di,vi,fc,fci(k),fs,B,Bi(k),NT,NTsw,Ti(k),Tswi(k),M,Mi(k),offset(k));
    y_tf_ts = y_target + yi + noise;
    sin_i(:,k) = (xi - conj(xi))/(1i*2);
    R(k) = max(xcorr(sin_t,sin_i(:,k),'normalized'),[],'all');
    if isnan(R(k))
        R(k) = 0;
    end
    Pi(k) = sum(yi.*conj(yi),"all");
    SIR(k) = Pt/(Pi(k) + Pn);
    SIRdB(k) = 10*log10(Pt/(Pi(k) + Pn));
    y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);
    v_ax = ((-M/2):(M/2))*delta_v;
    d_ax = (0:(N-1))*delta_d;

    if dB_map_flag ~= 0
        figure
        imagesc(v_ax(round(M/4):round(3*M/4)),...
            flipud(d_ax(1:round(size(d_ax,2)/2))'),...
            flipud(20*log10(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))/...
            max(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))),[],'all')))))
        xlabel('Velocity [m/s]','FontSize',30)
        ylabel('Distance [m]','FontSize',30)
        cb = colorbar;
        cb.Label.String = ('[dB]');
        cb.FontSize = 30;
        set(gca,'YDir','normal')
        title({['SIR = ', num2str(SIRdB(k)), ' dB | Correlation = ', num2str(R(k)) ] ...
            ['Offset = ', num2str(offset(k)) ,' | fci/fc = ', num2str(fci(k)/fc), ' | Bi/B = ', num2str(Bi(k)/B)] ...
            ['Up Ramp time ratio = ' , num2str(Ti(k)/(NT*Ts)), ' | \DeltaTi/\DeltaT = ' , num2str((Tswi(k)-Ti(k))/((NTsw-NT)*Ts))] ...
            ['N. of ramps in one CPI= ' , num2str(Mi(k))]}, 'FontSize', 30);
        set(gca,'fontsize',30)
    end

    if linear_map_flag ~= 0
        figure
        imagesc(v_ax(round(M/4):round(3*M/4)),...
            flipud(d_ax(1:round(size(d_ax,2)/2))'),...
            flipud(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))/...
            max(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))),[],'all'))))
        xlabel('Velocity [m/s]','FontSize',30)
        ylabel('Distance [m]','FontSize',30)
        cb = colorbar;
        cb.FontSize = 30;
        set(gca,'YDir','normal')
        title({['SIR = ', num2str(SIRdB(k)), ' dB | Correlation = ', num2str(R(k)) ] ...
            ['Offset = ', num2str(offset(k)) ,' | fci/fc = ', num2str(fci(k)/fc), ' | Bi/B = ', num2str(Bi(k)/B)] ...
            ['Up Ramp time ratio = ' , num2str(Ti(k)/(NT*Ts)), ' | \DeltaTi/\DeltaT = ' , num2str((Tswi(k)-Ti(k))/((NTsw-NT)*Ts))] ...
            ['N. of ramps in one CPI= ' , num2str(Mi(k))]}, 'FontSize', 30);
        set(gca,'fontsize',30)
    end
end
