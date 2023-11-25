clear all
clc
close 

%% Costanti
c = 3e8; % Velocità della luce

%% Bersagli
N_targets = 1;
d_targets = [150]; % Distanza dei bersagli in metri
v_targets = [50]; % Velocità dei bersagli [m/s]
v_targets_km_h = v_targets*3600/1000; % Velocità dei bersagli [km/h]
if (N_targets ~= size(d_targets,2) || N_targets ~= size(v_targets,2))
    disp('Errore nel fornire le caratteristiche dei bersagli. Controllare che velocità e distanze fornite siano congruenti al numero di bersagli.')
    return
end

%%  Parametri di sistema 
du = 600; % Massima distanza non ambigua
dmax = 600; % Massima distanza di interesse
vu_i = 115; % Massima velocità non ambigua ideale [m/s]
Nc = 512; % Numero di subcarrier
M = 4; % Ordine della modulazione da usare sulle sottoportanti
% fc = 140e9; % Frequenza della portante
fc = 140e9; % Frequenza della portante
pilots = [1 Nc/2 Nc]'; % Toni pilota

%% Parametri derivati
lambda_c = freq2wavelen(fc); % Lunghezza d'onda della portante
delta_f = c/(2*du); % Spaziatura fra le portanti
N_ofdm = Nc*4; % Numero di simboli OFDM da inviare in un ciclo di misurazione
B = delta_f*Nc; % Banda del segnale OFDM
T = 1/delta_f; % Durata simbolo OFDM
T_s = 1/B; % Periodo di campionamento del segnale
T_cp_i = 2*dmax/c; % Durata prefisso ciclico
N_cp = ceil(T_cp_i/T_s); % Campioni del prefisso ciclico
N_cp = 0;
T_cp = N_cp*T_s; % Durata effettiva del prefisso ciclico
T_ofdm = T_cp + T; % Durata del simbolo OFDM compreso di prefisso ciclico
T_sri_noT0 = N_ofdm*T_ofdm;
vu_noT0 = c /(2*fc*T_sri_noT0);
if vu_noT0 < vu_i
    T_sri_i = c/(2*fc*vu_i);
    N_0 = ceil(T_sri_i-T_sri_noT0); % Lunghezza dello zero-padding per raggiungere il requisito sulla massima velocità non ambigua
else
    N_0 = 0; % Non c'è bisogno di zero-padding perché i requisiti sono raggiunti
end
N_0 = Nc*4;
T_0 = N_0*T_s; % Intervallo di silenzio fra i simboli OFDM (durata dello zero-padding)
T_sri = T_ofdm + T_0; % Tempo di ripetizione del simbolo
T_cycle = N_ofdm*T_sri; % Durata di un ciclo di misurazione
vu = c/(2*fc*T_sri); % Massima velocità non-ambigua rilevabile [m/s]
delta_v = c/(2*fc*T_cycle); % Risoluzione in velocità
delta_fd = 2*delta_v*fc/c; % Risoluzione in frequenza Doppler
N_pilots = numel(pilots); % Numero di portanti pilota
N_dati_pilots = N_pilots*N_ofdm; % Numero di dati da inviari sulle portanti pilota
N_dati = (Nc-N_pilots)*N_ofdm; % Numeri di dati modulati da inviare
delta_d = c/(2*B); % Risoluzione in range

%% Generazione dati no piloti
% data = randi([0 M-1], Nc);
% parallel_mod_data = pskmod(data, M);
% serial_all_data = reshape(data,Nc*Nc,1);


%% Generazione dei dati
data = randi([0 M-1], N_dati, 1);
data_pilots = ones(N_pilots,N_ofdm);

%% Modulazione dei dati
mod_data = pskmod(data, M, pi/M);
mod_pilots = pskmod(data_pilots, M, pi/M);

%% Serial-to-Parallel e inserimento portanti pilota
parallel_mod_data = reshape(mod_data, Nc-N_pilots, []);
if N_pilots>0
    for i = 1:N_pilots
        if pilots(i) == 1
            parallel_mod_data = [mod_pilots(i,:) ; parallel_mod_data];
        elseif pilots(i) == Nc
            parallel_mod_data = [parallel_mod_data ; mod_pilots(i,:)];
        else
            parallel_mod_data = [parallel_mod_data(1:pilots(i)-1,:) ; mod_pilots(i,:) ; parallel_mod_data(pilots(i):end,:)];
        end
    end
end

serial_all_data = pskdemod(reshape(parallel_mod_data,Nc*N_ofdm,1),M);

%% Calcolo di x_t senza ifft, aggiunta cp e zero padding
t = 0:T_s:Nc*T;
t(end) = [];
t = t';
t_cp = 0:T_s:(Nc+N_cp)*T;
t_cp(end) = [];
t_cp = t_cp';
x_t = zeros(Nc, N_ofdm);

for u = 1:Nc
    for n = 1:Nc
        x_t(:,u) = x_t(:,u) + parallel_mod_data(n,u)*exp(1j*2*pi*(n-1)*delta_f*t((u-1)*Nc+1:u*Nc));
    end
end
x_t = x_t/Nc;
x_t_cp = [x_t(end-N_cp+1:end,:) ; x_t];
zp = zeros(N_0,N_ofdm);
x_t_cp_0 = [x_t_cp ; zp];
x_t_serial = reshape(x_t_cp_0, N_ofdm*(Nc+N_cp+N_0),1);

%% Quadrature modulation
t_tot = 0:T_s:T_cycle;
t_tot(end) = [];
t_tot = t_tot';
x_t_qam = x_t_serial./exp(1j*2*pi*fc*t_tot);

%% Generazione segnale ricevuto
L = zeros(N_targets,1);
for i = 1:N_targets
    L(i) = sqrt(10^(-fspl(d_targets(i),lambda_c)/20));
end
r_t = zeros(numel(t_tot),N_targets); % Distanza dei bersagli nel tempo
delay_t = zeros(numel(t_tot),N_targets); % Ritardo dei segnali riflesso dai bersagli nel tempo
for i = 1:N_targets
    r_t(:,i) = d_targets(i)+t_tot*v_targets(i);
    delay_t(:,i) = 2*r_t(:,i)/c;
end

echoes_m = zeros(Nc, N_ofdm, N_targets);
for l = 1:N_targets
    for u = 1:Nc
        for n = 1:Nc
            echoes_m(:,u,l) = echoes_m(:,u,l) + L(l)*parallel_mod_data(n,u)*exp(1j*2*pi*(n-1)*delta_f*(t_tot(u*N_cp+(u-1)*Nc+(u-1)*N_0+1:u*N_cp+u*Nc+(u-1)*N_0,1)-delay_t(u*N_cp+(u-1)*Nc+(u-1)*N_0+1:u*N_cp+u*Nc+(u-1)*N_0,l)));
        end
    end
end

echoes_m = [echoes_m(end-N_cp+1:end,:,:) ; echoes_m ; zeros(N_0,N_ofdm,N_targets)];
echoes_t = reshape(echoes_m,N_ofdm*(Nc+N_cp+N_0),1,N_targets);
for l = 1:N_targets
    echoes_t(:,:,l) = echoes_t(:,:,l).*exp(1j*2*pi*fc*(t_tot-delay_t(:,l)));
end

y_t_nonoise = sum(echoes_t,3);
noise = randn(size(y_t_nonoise,1),size(y_t_nonoise,2))/5e4;
noise = zeros(size(y_t_nonoise,1),size(y_t_nonoise,2));
y_t = y_t_nonoise + noise;
snr_rx = snr(y_t, noise);

%% Battimento, Serial-to-Parallel e rimozione CP
y_t = y_t.*exp(-1j*2*pi*fc*(t_tot));
y_tf_ts_cp_t0 = reshape(y_t,Nc+N_cp+N_0,N_ofdm);
y_tf_ts = y_tf_ts_cp_t0(N_cp+1:N_cp+Nc,:);

%% FFT 
y_f_ts = fft(y_tf_ts,Nc,1);

%% Divisione simboli
z_f_ts = y_f_ts./parallel_mod_data;

%% Elaborazione range 
z_f_v = fft(z_f_ts,N_ofdm,2);
z_f_v = fftshift(z_f_v,2);

%% Elaborazione Doppler
z_d_v = ifft(z_f_v,Nc,1);

v_ax = (-(N_ofdm)/2:(N_ofdm)/2)*delta_v;
d_ax = (0:Nc-1)*delta_d;

%% Mappa range-Doppler
[~,indx_dmax] = min(abs(d_ax-dmax));

imagesc(v_ax,flipud(d_ax(1:indx_dmax)'),fliplr(flipud(20*log10(abs(z_d_v((1:indx_dmax),:)).^2))))
% imagesc(v_ax,flipud(d_ax(1:indx_dmax)'),fliplr(flipud((abs(z_d_v((1:indx_dmax),:)).^2))))
xlabel('Velocity [m/s]')
ylabel('Distance [m]')
cb = colorbar;
% cb.Label.String = ('[dB]');
set(gca,'YDir','normal')
