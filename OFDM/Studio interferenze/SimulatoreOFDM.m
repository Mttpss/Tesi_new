clear all
close all

%% SETTING RESULTS WANTED

plotter_flag = 1;

%% SETTING RADAR PARAMETERS

loadRadarSetup_flag = 0;
if loadRadarSetup_flag ~= 0
    load('RadarParameters\FMCWParameterSetup_1.mat');
else
    du_ideal = 300;
    vu_ideal = 95;
    delta_d_ideal = 0.3/3.6;
    delta_v_ideal = 1/3.6;
    fc = 140e9;
    rs_flag = 1;
    
    [du,di,delta_d,vu,delta_v,delta_fd,delta_f,Nc,B,Ts,T,T0,N0,...
    Tcp,Ncp,Tsri,Nofdm,Tcycle,max_fd,max_fd_n,max_rm,max_rm_n,max_dm,...
    max_dm_n] = parametrization(rs_flag,du_id,di_id,delta_d_id,vu_id,...
    delta_v_id,fc);

end

c = 3e8;
lambda_c = c/fc;
fs = 1/Ts;
t = 0:Ts:(M*NTsw*Ts-Ts);
t_up = 0:Ts:((NT-1)*Ts);

%% TARGET'S ECHOES GENERATION

d = 100;
v = 0;
[y_target,x] = echoSingleTarget(d,v,fc,delta_f,Nc,Nofdm,T,Tsri,s);
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

%% INTERFERENCE POINT SETUP

load_point_flag = 0;

if load_point_flag == 1
    load('point.mat');
else
    numofpoint = 202; % must be even
    fci_flag = 0;
    Nci_flag = 0;
    Nofdmi_flag = 0;
    Tcpi_flag = 0;
    T0i_flag = 0;
    offseti_flag = 1;

    if fci_flag ~= 0
        fci = linspace(fc-B,fc,numofpoint/2);
        fci(end) = [];
        fci = [fci linspace(fc,fc+B,numofpoint/2)];
    else
        fci = ones(1,numofpoint-1)*fc;
    end

    if Bi_flag ~= 0
        Bi = linspace(0,B,numofpoint/2);
        Bi(end) = [];
        Bi = [Bi linspace(B,B*2,numofpoint/2)];
    else
        Bi = ones(1,numofpoint-1)*B;
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
di = 0;
vi = 0;

R = zeros(numofpoint-1,1);
SIRdB = zeros(numofpoint-1,1);

for k = 1:(numofpoint-1)
    % disp('Iterazione ')
    % disp(k)
    [yi,xi] = echoInterferenceOFDM(x,Pt,n,di,vi,fc,fci(k),fs,B,Bi(k),NT,NTsw,Ti(k),Tswi(k),M,Mi(k),offset(k));
    % [yi,xi] = echoInterferenceFMCW(x,n,di,vi,fc,fc,fs,B,B,NT,NTsw,NT*Ts,NTsw*Ts,M,M,offset(k));
    y_i_noise = yi + noise;
    y_tf_ts = y_target + yi + noise;
    sin_i = (xi - conj(xi))/(1i*2);
    R(k) = xcorr(sin_t,sin_i,'normalized',0);
    % R(k) = max(abs(xcorr(x,xi,'normalized')),[],'all');
    if isnan(R(k))
        R(k) = 0;
    end
    Pi(k) = sum(yi.*conj(yi),"all");
    SIRdB(k) = 10*log10(Pt/(Pi(k) + Pn));
end

%% PLOTTING VARIOUS RESULTS

if plotter_flag ~= 0
    if fci_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,fci/fc,R,'--', 'Marker', 'o')
        title(ax1,'Carrier frequency ratio - Correlation')
        ax1.FontSize = 30;

        ax2 = nexttile;
        plot(fci/fc,SIRdB,'--', 'Marker', 'o')
        title(ax2,'Carrier frequency ratio - SIR')
        ax2.FontSize = 30;
    end

    if Bi_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Bi/B,R,'--', 'Marker', 'o')
        title(ax1,'Band ratio - Correlation')
        ax1.FontSize = 30;

        ax2 = nexttile;
        plot(ax2,Bi/B,SIRdB,'--', 'Marker', 'o')
        title(ax2,'Band ratio - SIR')
        ax2.FontSize = 30;
    end

    if Mi_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Mi,R,'--', 'Marker', 'o')
        title(ax1,'N. of ramp - Correlation')
        ax1.FontSize = 30;

        ax2 = nexttile;
        plot(ax2,Mi,SIRdB,'--', 'Marker', 'o')
        title(ax2,'N. of ramp - SIR')
        ax2.FontSize = 30;
    end

    if Ti_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Ti/(NT*Ts),R,'--', 'Marker', 'o')
        title(ax1,'Time of up ramp ratio - Correlation')
        ax1.FontSize = 30;

        ax2 = nexttile;
        plot(ax2,Ti/(NT*Ts),SIRdB,'--', 'Marker', 'o')
        title(ax2,'Time of up ramp ratio - SIR')
        ax2.FontSize = 30;
    end

    if Tswi_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,(Tswi-Ti)/((NTsw-NT)*Ts),R,'--', 'Marker', 'o')
        title(ax1,'\DeltaT ratio - Correlation')
        ax1.FontSize = 30;

        ax2 = nexttile;
        plot(ax2,(Tswi-Ti)/((NTsw-NT)*Ts),SIRdB,'--', 'Marker', 'o')
        title(ax2,'\DeltaT ratio - SIR')
        ax2.FontSize = 30;
    end

    if offset_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,offset,R,'--', 'Marker', 'o')
        title(ax1,'Offset - Correlation')
        ax1.FontSize = 30;

        ax2 = nexttile;
        plot(ax2,offset,SIRdB,'--', 'Marker', 'o')
        title(ax2,'Offset - SIR')
        ax2.FontSize = 30;
    end

    figure
    plot(R(ceil(end/2):end),SIRdB(ceil(end/2):end),'--', 'Marker', 'o')
    title('Correlation - SIR')
    fontsize(30,"points")
    % p.FontSize = 30;
end




