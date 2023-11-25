clear all
fc = 140e9;
rs_flag = 1;
du_id = 300;
di_id = du_id;
delta_d_id = 0.3;
vu_id = 360/3.6;
delta_v_id = 1/3.6;

[du,di,delta_d,vu,delta_v,delta_fd,delta_f,Nc,B,T_s,T,T0,N0,...
    Tcp,Ncp,Tsri,Nofdm,Tcycle,max_fd,max_fd_n,max_rm,max_rm_n,max_dm,...
    max_dm_n] = parametrizationOFDM(rs_flag,du_id,di_id,delta_d_id,vu_id,...
    delta_v_id,fc);

d = [100];
v = [40];

M = 4;
data = ones(Nc, Nofdm);
s = pskmod(data, M);

y_tf_ts = echoTarget(d,v,fc,delta_f,Nc,Nofdm,T,Tsri,s);

w = eye(Nc);
z_d_v = RDP2dClassic(y_tf_ts,s,w,Nc,Nofdm);

v_ax = (-(Nofdm)/2:(Nofdm)/2)*delta_v;
d_ax = (0:Nc-1)*delta_d;

[~,indx_dmax] = min(abs(d_ax-di_id));

%imagesc(v_ax,flipud(d_ax'),flipud(20*log10((abs(z_d_v)./max(abs(z_d_v),[],"all"))))) % Normalized
imagesc(v_ax,flipud(d_ax'),flipud((abs(z_d_v)./max(abs(z_d_v),[],"all")))) % Normalized
% imagesc(v_ax,flipud(d_ax'),flipud((abs(z_d_v))))

xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')