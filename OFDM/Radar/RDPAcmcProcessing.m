function [z_d_v] = RDPAcmcAcdcProcessing(y_tf_v,s,w_Nc,Nc,Nofdm)
% RDPAcmcAcdcProcessing is the function to perform range and Doppler processing on 
% the incoming echo with only ACMC performed before

y_f_v = fft(y_tf_v,Nc,1);
z_f_v = y_f_v./(s*ones(1,Nofdm));
z_f_v = fftshift(z_f_v,2);
z_d_v = ifft(w_Nc*z_f_v,Nc,1);

end