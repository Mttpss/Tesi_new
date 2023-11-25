function [y_d_v] = rangeDopplerProcessing(y_tf_ts,N,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y_d_ts = fft(y_tf_ts,N,1)./N;
y_d_v = fft(y_d_ts,M,2)./M;
y_d_v = fftshift(y_d_v,2);
end