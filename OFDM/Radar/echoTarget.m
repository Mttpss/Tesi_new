function [y_tf_ts] = echoTarget(d,v,fc,delta_f,Nc,Nofdm,T,Tsri,s)
%ECHOTARGET Summary of this function goes here
%   Detailed explanation goes here
if size(d) == size(v)
    Ntarget = max(size(d));
else
    display('Size of d and v are not compatible.');
end

y_tf_ts = zeros(Nc,Nofdm);
for i = 1:Ntarget
    y_tf_ts = y_tf_ts + echoSingleTarget(d(i),v(i),fc,delta_f,Nc,Nofdm,T,Tsri,s);
end
end

