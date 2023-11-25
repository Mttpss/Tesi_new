function [y,x] = echoInterferenceFMCW(x_tobeat,n,di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset)
Ts = 1/fs;

y = zeros(NT,M);
x = zeros(NTsw*M,1);

for i = 1:n
    [yi,xi] = echoSingleInterferenceFMCW(x_tobeat,di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset);
    y = y + yi;
    x = xi;
end

end