clear all
close all

%%% FMCW
NFFT=512;
NCHIRP=256;
FS=40e6;
TS=1/FS;
t=[0:TS:(NFFT-1)*TS];

R1=10; %m
R2=200; %m
V1=50/3.6;% m/s
V2=200/3.6; %m/s
f0=79e9; %Hz
t0=1./f0;

c=3e8; %m/s

DF=154e6;    %Chirpband B
Tsw=15.8e-6;    %Chirp time

m=DF/Tsw;

%%% per fare spttro banda chirp
NFFTT=2^14;
fc=50e6 %100 MHz
fc=77e6 %100 MHz
fof=50e6 % MHz

%FSS=770e6;TSS=1/FSS;
TSS=Tsw/NFFTT;FSS=1/TSS;
t1=[0:1:(NFFTT-1)]*Tsw/NFFTT;

st=sin(2*pi*(fc+0.5*m*t1).*t1);
ST=fft(st,NFFTT)/NFFTT;
fST=[0:1:NFFTT-1]/(NFFTT*TSS);
figure;plot(t1,st);
figure;plot(fST,abs(ST));
ylabel('Signal Amplitude [n.u.]','fontsize',18)
set(gca,'fontsize',18,'xlim',[0 FSS/2])
grid on

%%%%
%%%% noise
TK=293;
BN=FS;
SNR=0 %dB di rapposrto segnale rumore voluto
pn=1.380658e-23.*TK.*BN; %questo dovrebbe essere il KTB in Watt 

pn=0.5*10^(-SNR/10); %questo vale imponendo come riferimenteo di potenza per il seganale 0.5 (V^2/2 con V=1)

% genera i campioni di rumore AWGN
%rng('shuffle')
randn('state',sum(100*clock));
nn=sqrt(pn).*randn(size(t));
%   out_n=filtro_BP(f0,BLP,n,t,ts);   
%    nn=filtro_BP(1+FS/2,FS,nn,t,TS/3);   


T=[0:Tsw:(NCHIRP-1)*Tsw];
fd1=2*V1*f0/c+m*2*(R1+V1.*T)/c;
fv1=2*f0*(R1+V1.*T)/c;
fd2=2*V2*f0/c+m*2*(R2+V2.*T)/c;
fv2=2*f0*(R2+V2.*T)/c;

for i=1:NCHIRP 
    %x=sin(2*pi*(fd1(i).*t+fv1(i)))+sin(2*pi*(fd2(i).*t+fv2(i))); 
    x0=sin(2*pi*(fd1(i).*t+fv1(i)))+sin(2*pi*(fd2(i).*t+fv2(i))); %senza AWGN
    x=x0+nn; %con AWGN
    x20=sin(2*pi*((fof+fd1(i)).*t+fv1(i)))+sin(2*pi*(fd2(i).*t+fv2(i))); %senza AWGN
    x2=x20+nn; %con AWGN
    %x=sin(2*pi*fd1(i).*t+fv1(i)); 
    Y0(i,:)=fft(x0,NFFT)./NFFT; 
    Y(i,:)=fft(x,NFFT)./NFFT; 
    Y2(i,:)=fft(x2,NFFT)./NFFT; 
end
figure;plot(t,x,t,x2)
f = FS*linspace(0,1,NFFT);
R0=c.*f/(2*m);
%figure;plot(R0,abs(Y(1,:)));
figure;plot(f,abs(Y(1,:)));
%set(gca,'yscale','log')
title('')
xlabel('frequency [Hz]','fontsize',18)
ylabel('amplitude [n.u.]','fontsize',18)
set(gca,'xlim',[0,FS/2],'fontsize',18)
figure;plot(R0,abs(Y(1,:)));

for i=1:NFFT
    ZZ0(i,:)=fft(Y0(:,i),NCHIRP)./NCHIRP;
    ZZ(i,:)=fft(Y(:,i),NCHIRP)./NCHIRP;
    ZZ2(i,:)=fft(Y2(:,i),NCHIRP)./NCHIRP;
    
end

%ZZ=fft2(Y);
FSV=1/Tsw;
fsv = FSV*linspace(0,1,NCHIRP);
VD=3.6*(c*fsv/(2*f0));% km/h

[RR,VV]=meshgrid(R0,VD);
n=0;
%ZZ=fftshift(ZZ);
%ZZ2=fftshift(ZZ2);
%RR=RR-max(R0)/2;
%VV=VV-max(VD)/2;

figure;surface(RR,VV,abs(ZZ'));shading flat;colormap(1-gray);
n=n+1;h(n)=xlabel(['Range [m] step:',num2str(R0(2)-R0(1)),' m']);
n=n+1;h(n)=ylabel(['Speed [km/h] step:',num2str(VD(2)-VD(1)),' km/h']);
set(gca,'fontsize',18,'xlim',[0 max(max(RR))],'ylim',[-max(VD)/2 max(VD)/2]);

figure;surf(RR(1:NCHIRP/2,1:NFFT/2),VV(1:NCHIRP/2,1:NFFT/2),abs(ZZ(1:NFFT/2,1:NCHIRP/2)'));
n=n+1;h(n)=xlabel(['Range [m] step:',num2str(R0(2)-R0(1)),' m']);
n=n+1;h(n)=ylabel(['Speed [km/h] step:',num2str(VD(2)-VD(1)),' km/h']);
n=n+1;h(n)=zlabel(['Amplitude [n.u.]'])
title('ZZ')
set(gca,'fontsize',18,'xlim',[0 max(max(RR))]/2,'ylim',[0 max(VD)/2])

figure;surf(RR(1:NCHIRP/2,1:NFFT/2),VV(1:NCHIRP/2,1:NFFT/2),abs(ZZ0(1:NFFT/2,1:NCHIRP/2)'));
n=n+1;h(n)=xlabel(['Range [m] step:',num2str(R0(2)-R0(1)),' m']);
n=n+1;h(n)=ylabel(['Speed [km/h] step:',num2str(VD(2)-VD(1)),' km/h']);
n=n+1;h(n)=zlabel(['Amplitude [n.u.]'])
title('ZZ0')
set(gca,'fontsize',18,'xlim',[0 max(max(RR))]/2,'ylim',[0 max(VD)/2])

figure;surf(RR(1:NCHIRP/2,1:NFFT/2),VV(1:NCHIRP/2,1:NFFT/2),abs(Y0(1:NCHIRP/2,1:NFFT/2)));
n=n+1;h(n)=xlabel(['Range [m] step:',num2str(R0(2)-R0(1)),' m']);
n=n+1;h(n)=ylabel(['Speed [km/h] step:',num2str(VD(2)-VD(1)),' km/h']);
n=n+1;h(n)=zlabel(['Amplitude [n.u.]']);
set(gca,'fontsize',18,'xlim',[0 max(max(RR))]/2,'ylim',[0 max(VD)/2])

figure;surf(RR(1:NCHIRP/2,1:NFFT/2),VV(1:NCHIRP/2,1:NFFT/2),abs(Y(1:NCHIRP/2,1:NFFT/2)));
n=n+1;h(n)=xlabel(['Range [m] step:',num2str(R0(2)-R0(1)),' m']);
n=n+1;h(n)=ylabel(['Speed [km/h] step:',num2str(VD(2)-VD(1)),' km/h']);
n=n+1;h(n)=zlabel(['Amplitude [n.u.]']);
set(gca,'fontsize',18,'xlim',[0 max(max(RR))]/2,'ylim',[0 max(VD)/2])

figure;surface(RR,VV,abs(ZZ2'));shading flat;colormap(1-gray);
n=n+1;h(n)=xlabel(['Range [m] step:',num2str(R0(2)-R0(1)),' m'])
n=n+1;h(n)=ylabel(['Speed [km/h] step:',num2str(VD(2)-VD(1)),' km/h'])
set(gca,'fontsize',18,'xlim',[0 max(max(RR))],'ylim',[-max(VD)/2 max(VD)/2])
set(h,'fontsize',18)

title('Interferenza')

