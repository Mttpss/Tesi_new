function [y,x] = echoInterferenceFMCW(n,di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset);
Ts = 1/fs;

y = zeros(NT,M);
x = zeros(max((floor(Tswi/Ts).*Mi),1),1);
% f = zeros(max((floor(Tswi/Ts).*Mi),1),1);

for i = 1:n
    [yi,xi] = echoSingleInterferenceFMCW(di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset);
    y = y + yi;
    x(1:floor(Tswi(n)*Mi(i)/Ts),n) = xi;
    % f(1:floor(Tswi(n)/T),n) = fi;
end

end