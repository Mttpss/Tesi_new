function [y_tf_ts] = echoSingleTarget(d0,v,fc,delta_f,Nc,Nofdm,T,Tsri,s)
%%%% FUNCTION TO GENERATE OFDM ECHO FROM A SINGLE MOVING TARGET %%%%%%%%%%
% 
% Input parameter:
% - d0 = initial range of the target
% - v = velocity of the target
% - fc = frequency of the signal carrier
% - delta_f = subcarrier's spacing
% - Nc = number of subcarriers
% - Nofdm = number of OFDM in a CPI
% - T = time duration of a OFDM symbol
% - Tsri = Symbol Repetition Interval
% - s = symbol matrix
% 
% Output parameter:
% - y_tf_ts = received signal matrix in fast time and slow time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8; % Speed of light
lambda_c = c/fc; % Wavelength of the carrier

tau = 2*d0/c; % Delay of the echo
tau_n = tau/T; % Normalized delay

L = sqrt(10^(-fspl(2*d0,lambda_c)/20)); % Free space path loss of the echo
a = L*exp(-j*2*pi*fc*tau);

fd = 2*v*fc/c; % Doppler shift of the echo
fd_n = fd/delta_f; % Normalized Doppler shift

alpha = Tsri/T;

gamma = 2*v/c;

D1 = D_N_maker(fd_n/Nc,Nc);
D2 = conj(D_N_maker(tau_n,Nc));
D3 = D_N_maker(fd_n*alpha,Nofdm);

F_Nc = F_N_matrix_maker(Nc);
F_Nc_inv = inv(F_Nc);

P = P_maker(alpha,gamma,Nc,Nofdm);

y_tf_ts = a*D1*F_Nc_inv*D2*(s.*P)*D3;
end