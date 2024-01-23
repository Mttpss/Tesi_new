clear all
close all
%% SETTING FLAG FOR DESIRED MAPS
linear_map_flag = 0;
dB_map_flag = 1;

%% LOADING TARGETS AND INTERFERENCES
Id = 336;
moment = 0;
% loadfile = ["InterferingPosition/interferingPosition_Id" num2str(Id) "_moment" num2str(0) ".mat"];
load(sprintf("InterferingInfo/interferingPosition_Id%d_moment%d.mat",Id,moment));

%% SETTING RADAR PARAMETERS
Pt = 100; 
loadRadarSetup_flag = 0;
if loadRadarSetup_flag ~= 0
    load('RadarParameters\FMCWParameterSetup_1.mat');
else
    du_ideal = 300;
    vu_ideal = 150;
    delta_d_ideal = 0.3/3.6;
    delta_v_ideal = 1/3.6;
    B = 500e6;
    fc = 140e9;
    fs = 1300e6;
    %fs = 1000e6;

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

%% REFERENCE TARGET'S ECHOES GENERATION
d_ref = [0 100];
v_ref = [0 5];
[y_ref_target,x_ref] = echoSingleTarget(d_ref,v_ref,Pt,fc,fs,B,NT,NTsw,M);
PRef = sum(y_ref_target.*conj(y_ref_target),"all");
yv_ref = rangeDopplerProcessing(y_ref_target,N,M);
yv_ref_pow = yv_ref.*conj(yv_ref);
peakRefPow = max(yv_ref_pow,[],"all");
SNRtarget = 10*log10(sum(y_ref_target.*conj(y_ref_target),"all")/PRef);

%% TARGET'S ECHOES GENERATION
d = interferingInfo(:,1:2);
v = interferingInfo(:,3:4);
% d = [d ; 50 0 ; 100 0 ; 200 0];
% v = [v ; -1 0 ; -1 0 ; -1 0];
% d = [(50:10:200)' zeros(size((50:10:200)'))];
% v = rand(size(d))*2;
[y_target,Pecho,x] = echoTarget(d,v,Pt,fc,fs,B,NT,NTsw,M);
SNREcho = 10*log10(Pecho/PRef);

%% NOISE GENERATION
desiredSNR = 0;
noise = awgn(y_target,0) - y_target;
a = sqrt((PRef/(sum(noise.*conj(noise),"all")))*10^(-desiredSNR/10));
noise = a*noise;
noise_fft = fft(noise)/fs;
Pn = sum(noise.*conj(noise),"all");
SNRdB = 10*log10(PRef/Pn);

%% INTERFERENCE POINT SETUP

% load_point_flag = 0;
% 
% if load_point_flag == 1
%     load('point.mat');
% else
%     numofpoint = 20; % must be even
%     fci_flag = 0;
%     Bi_flag = 0;
%     Mi_flag = 0;
%     Ti_flag = 0;
%     Tswi_flag = 0;
%     offset_flag = 1;
% 
%     if Bi_flag ~= 0
%         Bi = linspace(0,B,numofpoint/2);
%         Bi(end) = [];
%         Bi = [Bi linspace(B,B*2,numofpoint/2)];
%     else
%         Bi = ones(1,numofpoint-1)*B;
%     end
% 
%     if fci_flag ~= 0
%         fci = linspace(fc-B,fc,numofpoint/2);
%         fci(end) = [];
%         fci = [fci linspace(fc,fc+B,numofpoint/2)];
%     else
%         fci = ones(1,numofpoint-1)*fc;
%     end
% 
%     if Mi_flag ~= 0
%         Mi = floor(linspace(0,M,numofpoint/2));
%         Mi(end) = [];
%         Mi = [Mi floor(linspace(M,2*M,numofpoint/2))];
%     else
%         Mi = ones(1,numofpoint-1)*M;
%     end
% 
%     if Ti_flag ~= 0
%         Ti = linspace(NT*Ts/10,NT*Ts,numofpoint/2);
%         Ti(end) = [];
%         Ti = [Ti linspace(NT*Ts,NT*Ts*2,numofpoint/2)];
%     else
%         Ti = ones(1,numofpoint-1)*NT*Ts;
%     end
% 
%     if Tswi_flag ~= 0
%         Tswi = floor(linspace(0,(NTsw-NT),numofpoint/2))*Ts;
%         Tswi(end) = [];
%         Tswi = [Tswi floor(linspace((NTsw-NT),(NTsw-NT)*2,numofpoint/2))*Ts];
%         Tswi = Ti + Tswi;
%     else
%         Tswi = ones(1,numofpoint-1)*(NTsw-NT)*Ts+Ti;
%     end
% 
%     if offset_flag ~= 0
%         offset = floor(linspace(-floor(max(Tswi,[],"all")/Ts*max(Mi,[],"all")),0,numofpoint/2));
%         offset(end) = [];
%         offset = [offset floor(linspace(0,NTsw*M,numofpoint/2))];
%     else
%         offset = zeros(1,numofpoint-1);
%     end
% end

%% INTERFERENCE GENERATION
n = size(interferingInfo,1);
di = interferingInfo(:,1:2);
vi = interferingInfo(:,3:4);
% di = [di ; 50 0 ; 100 0 ; 200 0];
% vi = [vi ; -1 0 ; -1 0 ; -1 0];
% n = n + 3;
[yi,Pi,xi] = echoInterferenceFMCW(x_ref,Pt,n,di,vi,fc,fc*ones(1,n),fs,B,B*ones(1,n),NT,NTsw,NT*Ts*ones(1,n),NTsw*Ts*ones(1,n),M,M*ones(1,n),zeros(1,n));
SNREcho_i = 10*log10(Pi/PRef);
y_tf_ts = y_ref_target + y_target + yi + noise;
% y_tf_ts = y_target + yi + noise;
% y_tf_ts = y_target + yi;
y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);

v_ax = ((-M/2):(M/2))*delta_v;
d_ax = (0:(N-1))*delta_d;

if dB_map_flag ~= 0
    figure
    imagesc(v_ax(round(M/4):round(3*M/4)),...
        flipud(d_ax(1:round(size(d_ax,2)/2))'),...
        flipud(10*log10(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4)).*conj(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4)))/peakRefPow)))
    % imagesc(v_ax(round(M/4):round(3*M/4)),...
    %     flipud(d_ax(1:round(size(d_ax,2)/2))'),...
    %     flipud(20*log10(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))/...
    %     max(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))),[],'all')))))
    xlabel('Velocity [m/s]','FontSize',30)
    ylabel('Distance [m]','FontSize',30)
    cb = colorbar;
    cb.Label.String = ('[dB]');
    cb.FontSize = 30;
    set(gca,'YDir','normal')
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
    set(gca,'fontsize',30)
end

load OtherTargetPosition_336_0.mat
otherTargetPosition = unique(otherTargetPosition,'rows');

figure
scatter(otherTargetPosition(:,1),otherTargetPosition(:,2),"filled","MarkerFaceColor","k","MarkerEdgeAlpha",.5)
set(gca,'YDir','normal')
set(gca,'fontsize',30)
axis([min([d(:,1) ; 0])-50 max([d(:,1) ; 0])+50 min([d(:,2) ; 0])-10 max([d(:,2) ; 0])+10])

hold on
scatter(d(:,1),d(:,2),200,"filled","MarkerFaceColor","r")
for i = 1:size(d,1)
    text(d(i,1),d(i,2), ['  \leftarrow v = ' num2str(v(i,1)) 'm/s'],'FontSize',25)
end

scatter(100,0,200,"filled","MarkerFaceColor","g")
text(100,0, ['  \leftarrow Reference target (distance = 100 m, speed = 5 m/s)'],'FontSize',25)

scatter(0,0,200,"filled","MarkerFaceColor","b")
text(0,0, ['  \leftarrow Interfered sensor'],'FontSize',25)

xlabel('x [m]','FontSize',40)
ylabel('y [m]','FontSize',40)
hold off



figure
scatter(otherTargetPosition(:,1),otherTargetPosition(:,2),"filled","MarkerFaceColor","k")
axis([min(otherTargetPosition(:,1)) max(otherTargetPosition(:,1)) min(otherTargetPosition(:,2)) max(otherTargetPosition(:,2))])
set(gca,'YDir','normal')
set(gca,'fontsize',30)
hold on
scatter(100,0,200,"filled","MarkerFaceColor","g")

scatter(0,0,200,"filled","MarkerFaceColor","b")

scatter(d(:,1),d(:,2),200,"filled","MarkerFaceColor","r")

xlabel('x [m]','FontSize',40)
ylabel('y [m]','FontSize',40)
hold off