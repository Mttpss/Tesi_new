function [y,x] = echoSingleInterferenceFMCW(Pt,di,vi,si,fc,fci,fs,delta_fi,Nc,Nci,Nofdm,Nofdmi,NTcp,NTcpi,NT0,NT0i,offset)
%ECHOSINGLETARGETFMCW Summary of this function goes here
%   Detailed explanation goes here

c = 3e8;
lambda_ci = c/fci;
Ts = 1/fs;
Ti = 1/delta_fi;
NTi = floor(Ti/Ts);
Ntot = (NTcpi + NTi + NT0i);
Ttot = Ntot*Ts;
t = (0:Ts:Ttot);
R = d + v*t;
L = sqrt(10.^(-fspl(2*abs(R),lambda_ci)/20)); % Free space path loss of the echo
L = reshape(L,Ntot,Nofdm);

% disp('x')
x = zeros(NTswi,Mi);

fni = (0:Nci)*delta_fi + fci;
x = zeros(Nc,(NTi+NT0i));

for u = 1:Nofdm
    for n = 1:Nc
        
    end
end

if offset < 0 
    x = [x; zeros(abs(offset),1)];
    x = x(end-NTswi*Mi+1:end);
elseif offset > 0
    x = [zeros(abs(offset),1) ; x];
    x = x(1:NTswi*Mi);
end
if size(x,1)*size(x,2) < NTsw*M
    x = [x ; zeros(NTsw*M - size(x,1)*size(x,2),1)];
else
    x = x(1:NTsw*M,1);
end

% figure
% spectrogram(x(1:NTsw*3),100,80,1024,fs,'yaxis')

y_prefilter = zeros(NTswi,Mi);
for l = 1:Mi
    y_prefilter(1:NTi,l) = L(1:NTi,l).*exp(1i*2*pi*(- fci*d/c - fci*v*(l-1)*NTswi*Ts/c + (fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c)*t_up + (mi/2 - mi*v/c)*(t_up.^2)))/sqrt(2);
    y_prefilter(NTi+1:NTswi,l) = zeros(NTswi-NTi,1);
end
a = sqrt((Pt/(sum(y_prefilter.*conj(y_prefilter),"all"))));
% y_prefilter = a*y_prefilter;
% SNI_prefilter_db = 10*log10(Pt/(sum(y_prefilter.*conj(y_prefilter),"all")));
% disp(SNI_prefilter_db)

y_filtered = zeros(NTswi,Mi);
if ((fci - mi*d/c - fci*v/c - mi*v*(Mi-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(end)))>=fc || (fci - mi*d/c - fci*v/c)<=fc+B
    for l = 1:Mi
        if ((fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(end)))>=fc || ((fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(1)))<=fc+B
            for k = 1:NTi
                if ((fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(k)))<fc || ((fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(k)))>fc+B
                    y_filtered(k,l) = 0;
                else
                    y_filtered(k,l) = L(k,l).*exp(1i*2*pi*(- fci*d/c - fci*v*(l-1)*NTswi*Ts/c + (fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c)*t_up(k) + (mi/2 - mi*v/c)*(t_up(k)^2)))/sqrt(2);
                end
            end
        end
        y_filtered((NTi+1):end,l) = zeros(NTswi-NTi,1);
    end
end
y_filtered = a*y_filtered;
% SNI_filtered_db = 10*log10(Pt/(sum(y_filtered.*conj(y_filtered),"all")));
% disp(SNI_filtered_db)
y_filtered = reshape(y_filtered,NTswi*Mi,1);

if offset < 0 
    y_filtered = [y_filtered; zeros(abs(offset),1)];
    y_filtered = y_filtered(end-NTswi*Mi+1:end);
elseif offset > 0
    y_filtered = [zeros(abs(offset),1) ; y_filtered];
    y_filtered = y_filtered(1:NTswi*Mi);
end

if size(y_filtered,1)*size(y_filtered,2) < NTsw*M
    y_filtered = [y_filtered ; zeros(NTsw*M - size(y_filtered,1)*size(y_filtered,2),1)];
else
    y_filtered = y_filtered(1:NTsw*M,1);
end

% figure
% spectrogram(y_prebeat_filtered(1:NTsw*3),100,80,1024,fs,'yaxis')

y = x_tobeat.*conj(y_filtered);
y = reshape(y,NTsw,M);
y = y(1:NT,:);
end

