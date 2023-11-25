clear all
desiredSNR = 10;
fc = 1e5;
fs = 1e6;
Ts = 1/fs;
t = 0:Ts:1;
seno_tempo = sin(2*pi*fc*t);
espo_tempo = exp(1i*2*pi*fc*t)/sqrt(2);
seno_frequenza = fft(seno_tempo)/fs;

noise_pre = awgn(seno_tempo,0) - seno_tempo;
noise_pre_filtered = noise_pre;
% noise_pre_filtered = lowpass(noise_pre,fc/fs);
a = sqrt((sum(seno_tempo.*conj(seno_tempo),2)/(sum(noise_pre_filtered.*conj(noise_pre_filtered),2)))*10^(-desiredSNR/10));
noise = a*noise_pre_filtered;
noise_fft = fft(noise)/fs;


SNR_tempo_seno = 10*log10(sum(seno_tempo.*conj(seno_tempo),2)/sum(noise.*conj(noise),2));
SNR_frequenza_seno = 10*log10(sum(seno_frequenza.*conj(seno_frequenza),2)/sum(noise_fft.*conj(noise_fft),2));
SNR_tempo_esp = 10*log10(sum(espo_tempo.*conj(espo_tempo),2)/sum(noise.*conj(noise),2));