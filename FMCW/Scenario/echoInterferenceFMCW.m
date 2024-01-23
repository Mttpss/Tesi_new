function [y,Pecho,x] = echoInterferenceFMCW(x_tobeat,Pt,n,di,vi,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset)
Ts = 1/fs;

y = zeros(NT,M);
x = zeros(NTsw*M,1);
Pecho = zeros(n,1);

for i = 1:n
    [yi,xi] = echoSingleInterferenceFMCW(x_tobeat,Pt,di(i,:),vi(i,:),fc,fci(i),fs,B,Bi(i),NT,NTsw,Ti(i),Tswi(i),M,Mi(i),offset(i));
    Pecho(i) = sum(yi.*conj(yi),"all");
    y = y + yi;
    x = xi;
end

end