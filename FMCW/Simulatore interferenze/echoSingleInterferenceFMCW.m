function [y,x] = echoSingleInterferenceFMCW(x_tobeat,d,v,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset)
%ECHOSINGLETARGETFMCW Summary of this function goes here
%   Detailed explanation goes here

%% PROVA
c = 3e8;
lambda_ci = c/fci;

fs_o = fs;


Ts = 1/fs;
NTi = floor(Ti/Ts);
NTswi = floor(Tswi/Ts);
mi = Bi/(NTi*Ts);

t = (0:Ts:(Mi*NTswi*Ts-Ts))';
t_up = (0:Ts:((NTi-1)*Ts))';
R = d + v*t;
L = sqrt(10.^(-fspl(2*R,lambda_ci)/20)); % Free space path loss of the echo
L = reshape(L,NTswi,Mi);

x = zeros(NTswi,Mi);
if fci + Bi > fc && fci < fc + B
    for l = 1:Mi
        %disp(['a' num2str(l)])
        x(1:NTi,l) = exp(1i*2*pi*(fci*t_up + mi/2*(t_up.^2)))/sqrt(2);
        x((NTi+1):end,l) = zeros(NTswi-NTi,1);
        % for k = 1:NTi
        %     if (fci+mi*t_up(k))<fc || (fci+mi*t_up(k))>fc+B
        %         x(k,l) = 0;
        %     end
        % end
    end
end
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

figure
spectrogram(x(1:NTsw*3),100,80,1024,fs_o,'yaxis')

y_prebeat = zeros(NTswi,Mi);

for l = 1:Mi
    y_prebeat(1:NTi,l) = L(1:NTi,l).*exp(1i*2*pi*(- fci*d/c - fci*v*(l-1)*NTswi*Ts/c + (fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c)*t_up + (mi/2 - mi*v/c)*(t_up.^2)))/sqrt(2);
    y_prebeat((NTi+1):end,l) = zeros(NTswi-NTi,1);
    for k = 1:NTi
        if ((fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(k)))<fc || ((fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c) + 2*(mi/2 - mi*v/c)*(t_up(k)))>fc+B
            y_prebeat(k,l) = 0;
        end
    end
end

y_prebeat = reshape(y_prebeat,NTswi*Mi,1);
y_prebeat_filtered = y_prebeat;

if offset < 0 
    y_prebeat_filtered = [y_prebeat_filtered; zeros(abs(offset),1)];
    y_prebeat_filtered = y_prebeat_filtered(end-NTswi*Mi+1:end);
elseif offset > 0
    y_prebeat_filtered = [zeros(abs(offset),1) ; y_prebeat_filtered];
    y_prebeat_filtered = y_prebeat_filtered(1:NTswi*Mi);
end

if size(y_prebeat_filtered,1)*size(y_prebeat_filtered,2) < NTsw*M
    y_prebeat_filtered = [y_prebeat_filtered ; zeros(NTsw*M - size(y_prebeat_filtered,1)*size(y_prebeat_filtered,2),1)];
else
    y_prebeat_filtered = y_prebeat_filtered(1:NTsw*M,1);
end

figure
spectrogram(y_prebeat_filtered(1:NTsw*3),100,80,1024,fs_o,'yaxis')

y = x_tobeat.*conj(y_prebeat_filtered);
y = reshape(y,NTsw,M);
y = y(1:NT,:);



%% VERA FUNZIONE
% c = 3e8;
% lambda_ci = c/fci;
% 
% fs_o = fs;
% i = 1;
% if fci < fc && (fci + Bi) < (fc + B)    
%     while i*fs_o < 2*(fc + B -fci)
%         i = i+1;
%     end
% elseif fci < fc && (fci + Bi) >= (fc + B)
%       while i*fs_o < 2*(Bi)
%         i = i+1;
%       end
% elseif fci >= fc && (fci + Bi) < (fc + B)
%     while i*fs_o < 2*(B)
%         i = i+1;
%     end
% elseif fci >= fc && (fci + Bi) >= (fc + B)
%      while i*fs_o < 2*(fci - fc + Bi)
%         i = i+1;
%      end
% end
% fs = fs_o*i;
% 
% Ts = 1/fs;
% NTi = floor(Ti/Ts);
% NTswi = floor(Tswi/Ts);
% mi = Bi/(NTi*Ts);
% 
% t = (0:Ts:(Mi*NTswi*Ts-Ts))';
% t_up = (0:Ts:((NTi-1)*Ts))';
% R = d + v*t;
% L = sqrt(10.^(-fspl(2*R,lambda_ci)/20)); % Free space path loss of the echo
% L = reshape(L,NTswi,Mi);
% 
% x = zeros(NTswi,Mi);
% for l = 1:Mi
%     x(1:NTi,l) = exp(1i*2*pi*(fci*t_up + mi/2*(t_up.^2)))/sqrt(2);
%     x((NTi+1):end,l) = zeros(NTswi-NTi,1);
% end
% x = reshape(x,NTswi*Mi,1);
% 
% % if fci + Bi < fc || fci > fc + B
% %     x = zeros(size(x));
% %     disp('A')
% % elseif fci >= fc && fci + Bi <= fc + B
% %     x = x;
% %     disp('B')
% % elseif fci <= fc && fci + Bi <= fc + B
% %     x = x.*exp(-1i*pi*fci*t);
% %     x = highpass(x,fc-fci,fs);
% %     x = x.*exp(1i*pi*fci*t); 
% %     disp('C')
% % elseif fci >= fc && fci + Bi >= fc + B 
% %     x = x.*exp(-1i*pi*fc*t);
% %     x = lowpass(x,B,fs);
% %     x = x.*exp(1i*pi*fc*t);
% %     disp('D')
% % elseif fci <= fc && fci + Bi >= fc + B
% %     x = x.*exp(-1i*pi*fci*t);
% %     x = lowpass(x,fc - fci + B,fs);
% %     x = highpass(x,fc-fci,fs);
% %     x = x.*exp(1i*pi*fci*t);
% %     disp('E')
% % end
% 
% x = decimate(x,i);
% if offset < 0 
%     x = [x; zeros(abs(offset),1)];
%     x = x(end-NTswi*Mi/i+1:end);
% elseif offset > 0
%     x = [zeros(abs(offset),1) ; x];
%     x = x(1:NTswi*Mi/i);
% end
% if size(x,1)*size(x,2) < NTsw*M
%     x = [x ; zeros(NTsw*M - size(x,1)*size(x,2),1)];
% else
%     x = x(1:NTsw*M,1);
% end
% figure
% spectrogram(x(1:NTsw*3),100,80,1024,fs/i,'yaxis')
% 
% y_prebeat = zeros(NTswi,Mi);
% for l = 1:Mi
%     y_prebeat(1:NTi,l) = L(1:NTi,l).*exp(1i*2*pi*(- fci*d/c - fci*v*(l-1)*NTswi*Ts/c + (fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c)*t_up + (mi/2 - mi*v/c)*(t_up.^2)))/sqrt(2);
%     y_prebeat((NTi+1):end,l) = zeros(NTswi-NTi,1);
% end
% y_prebeat = reshape(y_prebeat,NTswi*Mi,1);
% y_prebeat_filtered = y_prebeat;
% % if fci + Bi < fc || fci > fc + B
% %     y_prebeat_filtered = zeros(size(x));
% %     disp('A')
% % elseif fci >= fc && fci + Bi <= fc + B
% %     y_prebeat_filtered = y_prebeat;
% %     disp('B')
% % elseif fci <= fc && fci + Bi <= fc + B
% %     y_prebeat_filtered = y_prebeat.*exp(-1i*pi*fci*t);
% %     y_prebeat_filtered = highpass(y_prebeat_filtered,fc-fci,fs);
% %     y_prebeat_filtered = y_prebeat_filtered.*exp(1i*pi*fci*t); 
% %     disp('C')
% % elseif fci >= fc && fci + Bi >= fc + B 
% %     y_prebeat_filtered = y_prebeat.*exp(-1i*pi*fc*t);
% %     y_prebeat_filtered = lowpass(y_prebeat_filtered,B,fs);
% %     y_prebeat_filtered = y_prebeat_filtered.*exp(1i*pi*fc*t);
% %     disp('D')
% % elseif fci <= fc && fci + Bi >= fc + B
% %     y_prebeat_filtered = y_prebeat.*exp(-1i*pi*fci*t);
% %     y_prebeat_filtered = lowpass(y_prebeat_filtered,fc - fci + B,fs);
% %     y_prebeat_filtered = highpass(y_prebeat_filtered,fc - fci,fs);
% %     y_prebat_filtered = y_prebeat_filtered.*exp(1i*pi*fci*t);
% %     disp('E')
% % end
% 
% y_prebeat_filtered = decimate(y_prebeat_filtered,i);
% 
% if offset < 0 
%     y_prebeat_filtered = [y_prebeat_filtered; zeros(abs(offset),1)];
%     y_prebeat_filtered = y_prebeat_filtered(end-NTswi*Mi/i+1:end);
% elseif offset > 0
%     y_prebeat_filtered = [zeros(abs(offset),1) ; y_prebeat_filtered];
%     y_prebeat_filtered = y_prebeat_filtered(1:NTswi*Mi/i);
% end
% 
% if size(y_prebeat_filtered,1)*size(y_prebeat_filtered,2) < NTsw*M
%     y_prebeat_filtered = [y_prebeat_filtered ; zeros(NTsw*M - size(y_prebeat_filtered,1)*size(y_prebeat_filtered,2),1)];
% else
%     y_prebeat_filtered = y_prebeat_filtered(1:NTsw*M,1);
% end
% 
% figure
% spectrogram(y_prebeat_filtered(1:NTsw*3),100,80,1024,fs/i,'yaxis')
% 
% y = x_tobeat.*conj(y_prebeat_filtered);
% y = reshape(y,NTsw,M);
% y = y(1:NT,:);

end

