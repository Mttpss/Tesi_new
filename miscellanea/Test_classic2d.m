clear all
fc = 140e9;
rs_flag = 1;
du_id = 300;
di_id = du_id;
delta_d_id = 0.3;
vu_id = 95;
delta_v_id = 1;

aaFlag = 0;

[du,di,delta_d,vu,delta_v,delta_fd,delta_f,Nc,B,T_s,T,T0,N0,...
    Tcp,Ncp,Tsri,Nofdm,Tcycle,max_fd,max_fd_n,max_rm,max_rm_n,max_dm,...
    max_dm_n] = parametrization(rs_flag,du_id,di_id,delta_d_id,vu_id,...
    delta_v_id,fc);

d = 200;
v = 30;

M = 4;
data = randi([0 M-1], Nc, Nofdm);
s = pskmod(data, M);

y_tf_ts = echoSingleTarget(d,v,fc,delta_f,Nc,Nofdm,T,Tsri,s);

y_f_ts = fft(y_tf_ts,Nc,1);

w = eye(Nc);
z_d_v = rangeDopplerProcessing(y_f_ts,s,w,Nc,Nofdm,aaFlag);

v_ax = (-(Nofdm)/2:(Nofdm)/2)*delta_v;
d_ax = (0:Nc-1)*delta_d;

[~,indx_dmax] = min(abs(d_ax-di_id));
imagesc(v_ax,flipud(d_ax(1:indx_dmax)'),flipud(20*log10(abs(z_d_v((1:indx_dmax),:)).^2)))

xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')