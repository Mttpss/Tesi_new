clear all
close all

%% SETTING RESULTS WANTED

plotter_flag = 1;

%% SETTING RADAR PARAMETERS

loadRadarSetup_flag = 1;
if loadRadarSetup_flag ~= 0
    load('RadarParameters\FMCWParameterSetup_1.mat');
else
    du_ideal = 300;
    vu_ideal = 95;
    delta_d_ideal = 0.3/3.6;
    delta_v_ideal = 1/3.6;
    B = 500e6;
    fc = 140e9;
    fs = 500e6;

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
    numofpoint = 102; % must be even
    fci_flag = 0;
    Bi_flag = 1;
    Mi_flag = 0;
    Ti_flag = 0;
    Tswi_flag = 0;
    offset_flag = 0;

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
        Ti = linspace(Ts,NT*Ts,numofpoint/2);
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
SIRdB = zeros(numofpoint-1,1);

for k = 1:(numofpoint-1)
    disp('Iterazione ')
    disp(k)
    [yi,xi] = echoInterferenceFMCW(x,n,di,vi,fc,fci(k),fs,B,Bi(k),NT,NTsw,Ti(k),Tswi(k),M,Mi(k),offset(k));
    % [yi,xi] = echoInterferenceFMCW(x,n,di,vi,fc,fc,fs,B,B,NT,NTsw,NT*Ts,NTsw*Ts,M,M,offset(k));
    y_i_noise = yi + noise;
    y_tf_ts = y_target + yi + noise;
    sin_i = (xi - conj(xi))/(1i*2);
    R(k) = max(xcorr(sin_t,sin_i,'normalized'),[],'all');
    % R(k) = max(abs(xcorr(x,xi,'normalized')),[],'all');
    if isnan(R(k))
        R(k) = 0;
    end
    Pi = sum(yi.*conj(yi),"all");
    SIRdB(k) = 10*log10(Pt/(Pi + Pn));
end

%% PLOTTING VARIOUS RESULTS

if plotter_flag ~= 0
    if fci_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,fci/fc,R)
        title(ax1,'Carrier frequency ratio - Correlation')
        ax1.FontSize = 14;
        ax1.XColor = 'blue';

        ax2 = nexttile;
        plot(fci/fc,SIRdB)
        title(ax2,'Carrier frequency ratio - SIR')
        ax2.FontSize = 14;
        ax2.XColor = 'blue';
    end

    if Bi_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Bi/B,R)
        title(ax1,'Band ratio - Correlation')
        ax1.FontSize = 14;
        ax1.XColor = 'blue';

        ax2 = nexttile;
        plot(ax2,Bi/B,SIRdB)
        title(ax2,'Band ratio - SIR')
        ax2.FontSize = 14;
        ax2.XColor = 'blue';
    end

    if Mi_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Mi,R)
        title(ax1,'N. of ramp - Correlation')
        ax1.FontSize = 14;
        ax1.XColor = 'blue';

        ax2 = nexttile;
        plot(ax2,Mi,SIRdB)
        title(ax2,'N. of ramp - SIR')
        ax2.FontSize = 14;
        ax2.XColor = 'blue';
    end

    if Ti_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Ti/(NT*Ts),R)
        title(ax1,'Time of up ramp ratio - Correlation')
        ax1.FontSize = 14;
        ax1.XColor = 'blue';

        ax2 = nexttile;
        plot(ax2,Ti/(NT*Ts),SIRdB)
        title(ax2,'Time of up ramp ratio - SIR')
        ax2.FontSize = 14;
        ax2.XColor = 'blue';
    end

    if Tswi_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,Tswi,R)
        title(ax1,'Time beetween foolowing ramps - Correlation')
        ax1.FontSize = 14;
        ax1.XColor = 'blue';

        ax2 = nexttile;
        plot(ax2,Tswi,SIRdB)
        title(ax2,'Time beetween foolowing ramps - SIR')
        ax2.FontSize = 14;
        ax2.XColor = 'blue';
    end

    if offset_flag ~= 0
        figure
        tiledlayout(2,1)

        ax1 = nexttile;
        plot(ax1,offset,R)
        title(ax1,'Offset - Correlation')
        ax1.FontSize = 14;
        ax1.XColor = 'blue';

        ax2 = nexttile;
        plot(ax2,offset,SIRdB)
        title(ax2,'Offset - SIR')
        ax2.FontSize = 14;
        ax2.XColor = 'blue';
    end

    figure
    plot(R(ceil(end/2):end),SIRdB(ceil(end/2):end))
    title('Correlation - SIR')
    p.FontSize = 14;
    p.XColor = 'blue';
end




