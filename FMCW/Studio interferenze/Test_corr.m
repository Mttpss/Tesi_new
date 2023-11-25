clear all
close all
fs = 1e6;
Ts = 1/fs;
fc1 = 1e5;

t = 0:Ts:0.05;
a = exp(1i*2*pi*fc1*t);
b = sin(2*pi*fc1*t);
c = (a - conj(a))/(2*1i);
Ra = max(xcorr(a),[],'all');
Rb = max(xcorr(b),[],'all');
Rc = max(xcorr(c),[],'all');
[Ra_corr,lag] = xcorr(a);
Ra_corr_real = real(Ra_corr);
Rb_corr = xcorr(b);
Rc_corr = xcorr(c);
Ra_real = max(Ra_corr_real,[],'all');

% figure
% tiledlayout(4,1)
% 
% % Top plot
% ax1 = nexttile;
% plot(ax1,lag,Ra_corr)
% title(ax1,'Ra_corr')
% ax1.FontSize = 14;
% ax1.XColor = 'red';
% 
% % Bottom plot
% ax2 = nexttile;
% plot(ax2,lag,Ra_corr_real)
% title(ax2,'Ra_corr_real')
% ax2.FontSize = 14;
% ax2.XColor = 'blue';
% 
% % Bottom plot
% ax3 = nexttile;
% plot(ax3,lag,Rb_corr)
% title(ax3,'Rb_corr')
% ax3.FontSize = 14;
% ax3.XColor = 'yellow';
% 
% % Bottom plot
% ax4 = nexttile;
% plot(ax4,lag,Rc_corr)
% title(ax4,'Rc_corr')
% ax4.FontSize = 14;
% ax4.XColor = 'green';
% 
