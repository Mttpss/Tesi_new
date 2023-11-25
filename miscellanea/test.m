fc = 10;
phi0 = 0;
theta0 = 0;
T = 4;
B = 10;

M = 1;
Tsw = 20;
Ts = 1/((fc+B)*100);
NT = ceil(T/Ts);
NTsw = ceil(Tsw/Ts);

t = 0:Ts:(M*NTsw*Ts-Ts);
tMatrix = reshape(t,NTsw,M);
tMatrixNodown = tMatrix(1:NT,M);
x = zeros(NTsw,M);
y_tf_ts = zeros(NT,M);

m = B/T;
m_down = -B/(Tsw-T);
% for l = 1:M
%     y_tf_ts(:,l) = exp(1i*2*pi*fc*(tMatrixNodown(:,l)-2*R(:,l)/c)+...
%                    pi*m*(t(:,l).^2)+phi0-theta0);
% end
% for l = 1:M
%     x(1:NT,l) = exp(1i*2*pi*fc*(tMatrix(1:NT,l)+pi*m*(tMatrix(1:NT,l).^2)));
%     x((1+NT):end,l) = exp(1i*2*pi*fc*(tMatrix((1+NT):end,l)+pi*m_down*(tMatrix((1+NT):end,l).^2)));
% end

for l = 1:M
    x(1:NT,l) = sin(2*pi*(tMatrix(1:NT,l)+pi*m*(tMatrix(1:NT,l).^2)));
    x((1+NT):end,l) = sin(2*pi*(tMatrix((1+NT):end,l)+pi*m_down*(tMatrix((1+NT):end,l).^2)));
end

y_linear = reshape(y_tf_ts,NT*M,1);
x_linear = reshape(x,NTsw*M,1);
x_phase = angle(x_linear);
x_freq = diff(x_phase)/(Ts*2*pi);
plot(t(1:(end-1)),x_freq)
