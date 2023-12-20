clear all


d = [100 250];
v = [30 0];

du_i = 300;
vu_i = 95;
delta_d_i = 0.3/3.6;
delta_v_i = 1/3.6;
B = 500e6;
fc = 140e9;
fs = 500e6;

[phyReal_flag,du,vu,delta_d,delta_v,Ts,Tsw,NTsw,T,NT,N,M] = ...
    parametrizationFMCW(du_i,vu_i,delta_d_i,delta_v_i,B,fc,fs);

if phyReal_flag == 0
    disp(['Non è possibile raggiungere contemporaneamente i requisiti' ...
        ' di massima distanza e velocità non-ambigue']);
    return;
end

y_tf_ts = echoTarget(d,v,fc,fs,B,NT,NTsw,M);
y_d_v = rangeDopplerProcessing(y_tf_ts,N,M);

v_ax = ((-M/2):(M/2))*delta_v;
d_ax = (0:(N-1))*delta_d;

imagesc(v_ax(round(M/4):round(3*M/4)),...
    flipud(d_ax(1:round(size(d_ax,2)/2))'),...
    flipud(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))/...
    max(abs(y_d_v((1:round(size(d_ax,2)/2)),round(M/4):round(3*M/4))),[],'all'))))



xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')