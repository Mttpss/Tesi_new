function [du,di,delta_d,vu,delta_v,delta_fd,delta_f,Nc,B,T_s,T,T0,N0,...
    Tcp,Ncp,Tsri,Nofdm,Tcycle,max_fd,max_fd_n,max_rm,max_rm_n,max_dm,...
    max_dm_n] = parametrization(rs_flag,du_id,di_id,delta_d_id,vu_id,...
    delta_v_id,fc)

%%%% FUNCTION FOR PARAMETRIZATION OF AN CP-OFDM OR A RS-OFDM RADAR %%%%%%%
% 
% Input parameters:
% - rs_flag = flag for Ciclic Prefix (0) or Repeated Symbol (1) OFDM radar
% - du_id = ideal max umabiguous range
% - di_id = ideal max range of interest
% - delta_d_id = ideal range resolution
% - vu_id = ideal max unambiguos velocity
% - delta_v_id = ideal velocity resolution
% - fc = frequency of the carrier
% 
% Output parameters:
% - du = actual max unambiguos range
% - di = actual max range of interest
% - delta_d = actual resolution in range
% - vu = actual max unambiguous velocity
% - delta_v = actual velocity resolution
% - delta_fd = Doppler resolution
% - delta_f = subcarrier spacing
% - Nc = number of subcarriers
% - B = band of the OFDM signal
% - T_s = sample frequency
% - T = time duration of one OFDM symbol
% - T0 = time interval between one OFDM symbol and the next
% - N0 = number of 0 sample between one OFDM symbol and the next
% - Tcp = time duration of the CP
% - Ncp = number of sample of the CP
% - Tsri = Symbol Repetition Interval
% - Nofdm = number of OFDM symbol in a Coherent Processing Interval (CPI)
% - Tcycle = time duration of a CPI
% - max_fd = max Doppler shift obtainable with such vu
% - max_fd_n = normalized max_fd with regards to delta_f
% - max_rm = max range migration obtainable with such vu
% - max_rm_n = normalized max_rm with regards to delta_d
% - max_dm = max Doppler migration obtainable with such vu
% - max_dm_n = normalized max_dm with regards to delta_fd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8; % Speed of light [m/s]

delta_f = c/(2*du_id); % Subcarrier spacing [Hz]
du = du_id;
Nc = 2^(nextpow2(c/(2*delta_d_id*delta_f))); % Number of subcarriers
B = delta_f*Nc; % OFDM signal's band [Hz]
T_s = 1/(2*B); % Sample frequency
T = 1/delta_f; % Single OFDM symbol duration
delta_d = c/(2*Nc*delta_f); % Effective range resolution
if (di_id == du_id)
    Tcp = 0;
    di = du;
else
    Tcp = 2*di_id/c;
    Ncp = ceil(Tcp/T_s); % CP's number of samples
    Tcp = Ncp*T_s; % CP's duration
    di = Tcp*c/2; % Effective maximum distance of interest
end
Ncp = ceil(Tcp/T_s); % CP's number of samples
Tcp = Ncp*T_s; % CP's duration

if (rs_flag == 1)
    T0 = c/(2*fc*vu_id) - T_s;
else
    T0 = c/(2*fc*vu_id) - T_s - Tcp;
end
N0 = ceil(T0/T_s); % Number of 0 sample between OFDM symbols
T0 = N0*T_s; % Wait time between OFDM symbols [s]
if (rs_flag == 1) % Duration of the Symbol Repetion Interval
    Tsri = T + T0;
else
    Tsri = Tcp + T + T0;
end
vu = c/(2*fc*Tsri); % Effective maximum unambiguous velocity

Nofdm = 2^nextpow2(c/(2*fc*Tsri*delta_v_id)); % Number of OFDM symbol of a Coherent Processing Interval
Tcycle = Nofdm*Tsri; % Duration of a CPI
delta_v = c/(2*fc*Tcycle); % Effective velocity resolution
delta_fd = 2*delta_v*fc/c; % Doppler resolution

max_fd = 2*vu*fc/c; % Max Doppler shift
max_fd_n = max_fd/delta_f; % Mad normalized Doppler shift

max_rm = vu*Tcycle; % Max range migration
max_rm_n = max_rm/delta_d; % Max normalized range migration

max_dm = 2*vu*B/c; % Max Doppler migration
max_dm_n = max_dm/delta_fd; % Max normalized Doppler migration

end