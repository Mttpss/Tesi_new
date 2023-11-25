function [z_d_v] = RDP2dClassic(y_tf_ts,s,w,Nc,Nofdm)
% RDP2dClassic is the function to perform range and Doppler processing on 
% the incoming echo without any correction step (such as ACMC or ACDM)
% performed before

y_f_ts = fft(y_tf_ts,Nc,1);
z_f_ts = y_f_ts./s;
z_f_v = fft(z_f_ts,Nofdm,2);
z_f_v = fftshift(z_f_v,2);
z_d_v = ifft(z_f_v,Nc,1);

end