function [y_f_v] = ACMC(y_f_ts,w,delta_v,alpha,Nc,Nofdm)
%ACMC Summary of this function goes here
% y_f_ts has to be provided already corrected from the complex amplitude:
%           y_tf_ts_noamp = y_ts_tf*(diag(conj(s_s)/(norm(s_s)^2)));
%           y_f_ts = fft(y_tf_ts_noamp,Nc,1);
c = 3e8;
y_f_v = zeros(Nc,Nofdm);
F_Nofdm = F_N_matrix_maker(Nofdm);
M_n = zeros(Nofdm);
for n = 0:(Nc-1)
    for l = 0:(Nofdm-1)
        for u = 0:(Nofdm-1)
            M_n(l+1,u+1) = exp(1i*2*pi*n*(2*l*delta_v/c)*u*alpha);
        end
    end
    y_f_v(n+1,:) = (y_f_ts(n+1,:)*w)*(M_n.*F_Nofdm);
end
end