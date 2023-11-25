function [z_d_v] = rangeDopplerProcessing(y,s,w,Nc,Nofdm,aaFlag)
%rangeDopplerProcessing Summary of this function goes here
%   Detailed explanation goes here
if aaFlag == 0
    y_f_ts = y;
    z_f_ts = y_f_ts./s;
    z_f_v = fft(z_f_ts,Nofdm,2);
    z_f_v = fftshift(z_f_v,2);
    z_d_v = ifft(z_f_v,Nc,1);
elseif aaFlag == 1
    y_tf_v = y;
    y_f_v = fft(y_tf_v,Nc,1);
    z_f_v = y_f_v./(s*ones(1,Nofdm));
    z_f_v = fftshift(z_f_v,2);
    z_d_v = ifft(w*z_f_v,Nc,1);
end
end

