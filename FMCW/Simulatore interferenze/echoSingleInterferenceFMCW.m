function [y,x] = echoSingleInterferenceFMCW(x_tobeat,Pt,d,v,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset)
%ECHOSINGLETARGETFMCW Summary of this function goes here
%   Detailed explanation goes here

c = 3e8;
lambda_ci = c/fci;
Ts = 1/fs;
NTi = floor(Ti/Ts);
NTswi = floor(Tswi/Ts);
mi = Bi/(NTi*Ts);

t = (0:Ts:(Mi*NTswi*Ts-Ts))';
t_up = (0:Ts:((NTi-1)*Ts))';
R = d + v*t;
L = sqrt(10.^(-fspl(2*abs(R),lambda_ci)/20)); % Free space path loss of the echo
L = reshape(L,NTswi,Mi);

% disp('x')
x = zeros(NTswi,Mi);

if fci + Bi > fc && fci < fc + B
    for k = 1:NTi
        if (fci+mi*t_up(k))<fc || (fci+mi*t_up(k))>fc+B
            x(k,1) = 0;
        else
            x(k,1) = exp(1i*2*pi*(fci*t_up(k) + mi/2*(t_up(k)^2)))/sqrt(2);
        end
    end
    x((NTi+1):end,1) = zeros(NTswi-NTi,1);
    for l = 2:Mi
        x(:,l) = x(:,1);
    end
end
% for l = 1:Mi
%     x(1:NTi,l) = exp(1i*2*pi*(fci*t_up + mi/2*(t_up.^2)))/sqrt(2);
%     x(NTi+1:NTswi,l) = zeros(NTswi-NTi,1);
% end

x = reshape(x,NTswi*Mi,1);

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

