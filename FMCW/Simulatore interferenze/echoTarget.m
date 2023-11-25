function [y_tf_ts] = echoTarget(d,v,fc,fs,B,NT,NTsw,M)
%ECHOTARGET Summary of this function goes here
%   Detailed explanation goes here
if size(d) == size(v)
    Ntarget = max(size(d));
else
    display('Size of d and v are not compatible.');
end

y_tf_ts = zeros(NT,M);

for i = 1:Ntarget
    y_tf_ts = y_tf_ts + echoSingleTarget(d(i),v(i),fc,fs,B,NT,NTsw,M);
end
end

