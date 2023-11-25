clear all
close all
load("x_o.mat")
fc = 1.4e11;
fci = 1.4e11-170e6;
B = 500e6;
Bi = 833e6;
fs = 2*B;
Ts = 1/fs;
NT = 2000;
NTsw = 2820;
Ti = NT*Ts;
Tswi = NTsw*Ts;
M = 1024;
Mi = 1024;
offset = 0;

fs_o = fs;
i = 1;
while i*fs_o < 2*Bi
    i = i+1;
end
fs = fs_o*i;

c = 3e8;
lambda_ci = c/fci;
Ts = 1/fs;
m = B/(NT*Ts);
NTi = floor(Ti/Ts);
NTswi = floor(Tswi/Ts);
mi = Bi/(NTi*Ts);

t = 0:Ts:(Mi*NTswi*Ts-Ts);
t_up = (0:Ts:((NTi-1)*Ts))';
t_linear = t';

t = reshape(t,NTswi,Mi);

x = zeros(NTswi,Mi);
f = zeros(NTswi,Mi);
for l = 1:Mi
    x(1:NTi,l) = exp(1i*2*pi*(fci*t_up + mi/2*(t_up.^2)))/sqrt(2);
    x((NTi+1):end,l) = zeros(NTswi-NTi,1);
    f(1:NTi,l) = fci + mi/2*t_up;
    f((NTi+1):end,l) = zeros(NTswi-NTi,1);
end
x = reshape(x,NTswi*Mi,1);
disp('filtering')
if fci + Bi < fc || fci > fc + B
    x_filtered = zeros(size(x));
elseif fci >= fc && fci + Bi <= fc + B
    x_filtered = x;
elseif fci <= fc && fci + Bi <= fc + B
    x_filtered = x.*exp(-1i*pi*fci*t_linear);
    x_filtered = highpass(x_filtered,fc-fci,fs);
    x_filtered = x_filtered.*exp(1i*pi*fci*t_linear); 
elseif fci >= fc && fci + Bi >= fc + B 
    x_filtered = x.*exp(-1i*pi*fc*t_linear);
    x_filtered = lowpass(x_filtered,B,fs);
    x_filtered = x_filtered.*exp(1i*pi*fc*t_linear);
elseif fci <= fc && fci + Bi >= fc + B
    x_filtered = x.*exp(-1i*pi*fci*t_linear);
    x_filtered = lowpass(x_filtered,B,fs);
    x_filtered = highpass(x_filtered,fc-fci,fs);
    x_filtered = x_filtered.*exp(1i*pi*fci*t_linear);
end

disp('resampling')
x = decimate(x,i);
x_filtered = decimate(x_filtered,i);
disp('offsetting')
if offset < 0 
    x = [x(abs(offset)+1:end) ; zeros(abs(offset),1)];
    x_filtered = [x_filtered(abs(offset)+1:end) ; zeros(abs(offset),1)];
elseif offset > 0
    x = [zeros(abs(offset),1) ; x(1:(NTsw*M-abs(offset)))];
    x_filtered = [zeros(abs(offset),1) ; x_filtered(1:(NTsw*M-abs(offset)))];
end

x = x(1:NTsw*M,1);
x_filtered = x_filtered(1:NTsw*M,1);
disp('plotting')
figure
spectrogram(x(1:NTsw*3),100,80,2048,fs_o,'yaxis')

figure
spectrogram(x_filtered(1:NTsw*3),100,80,2048,fs_o,'yaxis')

figure
spectrogram(x_o(1:NTsw*3),100,80,2048,fs_o,'yaxis')

figure
spectrogram(x(1:NTsw*3)+x_o(1:NTsw*3)+x_filtered(1:NTsw*3),100,80,2048,fs_o,'yaxis')



