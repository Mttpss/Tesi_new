function [y,x] = echoInterferenceFMCW(Pt,di,vi,fc,fci,fs,Nc,Nci,Nofdm,Nofdmi,NTcp,NTcpi,NT0,NT0i,offset);
Ts = 1/fs;

y = zeros(Nc/Ts,Nofdm);
x = zeros(Nc*Nofdm/Ts,1);

for i = 1:n
    [yi,xi] = echoSingleInterferenceOFDM(Pt,di,vi,fc,fci,fs,Nc,Nci,Nofdm,Nofdmi,NTcp,NTcpi,NT0,NT0i,offset);
    y = y + yi;
    x = xi;
end

end