clear all
c = 3e8;
fc = 140e9;
rs_flag = 1;
du_id = 300;
di_id = du_id;
delta_d_id = 0.3;
vu_id = 360/3.6; % Meglio inserire un valore leggermente superiore al requisito richiesto
delta_v_id = 1/3.6;

[du,di,delta_d,vu,delta_v,delta_fd,delta_f,Nc,B,T_s,T,T0,N0,...
    Tcp,Ncp,Tsri,Nofdm,Tcycle,max_fd,max_fd_n,max_rm,max_rm_n,max_dm,...
    max_dm_n] = parametrizationOFDM(rs_flag,du_id,di_id,delta_d_id,vu_id,...
    delta_v_id,fc);

d = [100];
v = [40];

M = 4;
data = randi([0 M-1], Nc, 1);
s_c = pskmod(data, M);
s_s = ones(Nofdm,1);
s = s_c*s_s';
w_Nc = eye(Nc);
w_Nofdm = eye(Nofdm);

y_tf_ts = echoTarget(d,v,fc,delta_f,Nc,Nofdm,T,Tsri,s);
y_f_ts = fft(y_tf_ts,Nc,1);
y_f_v = ACMC(y_f_ts,w_Nofdm,delta_v,Tsri/T,Nc,Nofdm);
y_tf_v = ifft(y_f_v,Nc,1);
y_tf_v = ACDC(y_tf_v,w_Nc,Nc,Nofdm,Tsri/T);

z_d_v = RDPAcmcAcdcProcessing(y_tf_v,s_c,w_Nc,Nc,Nofdm);

v_ax = (-(Nofdm)/2:(Nofdm)/2)*delta_v;
d_ax = (0:Nc-1)*delta_d;

% imagesc(v_ax,flipud(d_ax'),flipud(20*log10((abs(z_d_v)./max(abs(z_d_v),[],"all"))))) % Normalized
imagesc(v_ax,flipud(d_ax'),flipud((abs(z_d_v)./max(abs(z_d_v),[],"all"))))
% imagesc(v_ax,flipud(d_ax'),flipud((abs(z_d_v))))

xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')