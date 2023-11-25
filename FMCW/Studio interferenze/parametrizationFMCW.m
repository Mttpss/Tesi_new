function [phyReal_flag,du,vu,delta_d,delta_v,Ts,Tsw,NTsw,T,NT,N,M] = parametrizationFMCW(du_i,vu_i,delta_d_i,delta_v_i,B,fc,fs)
% PARAMETRIZATIONFMCW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
% - du_i = ideal max non-ambiguous range
% - vu_i = ideal max non-ambiguous velocity
% - delta_d_i = ideal range resolution
% - delta_v_i = ideal velocity resolution
% - B = band of the FMCW signal
% - fc = frequency of the carrier
% - fs = sampling frequency
% 
% Output parameters:
% - phyReal_flag = physical reachability of both du_i and vu_i
% - du = actual max non-ambiguous range
% - vu = actual max non-ambiguous velocity
% - delta_d = actual
% - delta_v = actual
% - Ts = sampling periodo
% - Tsw = time beetwen one chirp's start and the following
% - T = time for the up-ramp to cover all the signal's band
% - NTsw = number of sample in one Tsw
% - NTs = number of sample in one Ts
% - N = point of the fft along the fast time axis
% - M = point of the fft along the slow time axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 3e8;

Ts = 1/fs;
Tsw = c/(4*fc*vu_i);
T = du_i*4*B*Ts/c;

if c*Tsw/2 < c*T/(4*B*Ts)
    phyReal_flag = 0;
    return;
else
    phyReal_flag = 1;
end

NTsw = ceil(Tsw/Ts);
vu = c/(4*fc*NTsw*Ts);

NT = round(T/Ts);
du = NT*Ts*c/(4*B*Ts);

Ntemp = c*T/(2*B*delta_d_i*Ts);
N = 2^nextpow2(Ntemp);
delta_d = c*T/(2*B*N*Ts);

Mtemp = c/(2*fc*delta_v_i*Tsw);
M = 2^nextpow2(Mtemp);
delta_v = c/(2*fc*M*Tsw);

end

