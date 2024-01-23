function [y_tf_ts,Pecho,x] = echoTarget(d,v,Pt,fc,fs,B,NT,NTsw,M)
%ECHOTARGET Summary of this function goes here
%   Detailed explanation goes here
if size(d,1) == size(v,1)
    NTarget = max(size(d));
else
    display('Size of d and v are not compatible.');
end

y_tf_ts = zeros(NT,M);
Pecho = zeros(NTarget,1);

for i = 1:NTarget
    [y,x] = echoSingleTarget(d(i,:),v(i,:),Pt,fc,fs,B,NT,NTsw,M);
    Pecho(i) = sum(y.*conj(y),"all");
    y_tf_ts = y_tf_ts + y;
end
end

