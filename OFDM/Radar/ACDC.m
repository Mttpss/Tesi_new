function [y_tf_v_corrected] = ACDC(y_tf_v_noamp,w,Nc,Nofdm,alpha)
%ACDC Summary of this function goes here
% y_tf_v has to be provided already corrected from the complex amplitude:
%           y_tf_ts_noamp = y_ts_tf*(diag(conj(s_s)/(norm(s_s)^2)));
%           y_tf_v = fft(y_tf_ts_noamp*w,Nofdm,2);
C = zeros(Nc,Nofdm);
for m = 0:(Nc-1)
    for l = 0:(Nofdm-1)
        C(m+1,l+1) = exp(-1i*2*pi*(l/(alpha*Nofdm))*m/Nc);
    end
end
y_tf_v_corrected = C.*y_tf_v_noamp;

end

