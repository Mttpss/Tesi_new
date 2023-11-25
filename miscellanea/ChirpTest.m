clear all
fc = 100;
phi0 = 0;
theta0 = 0;
T = 4;
B = 10;

M = 2;
Tsw = 50;
Ts = 1/((fc+B)*100);
NT = ceil(T/Ts);
NTsw = ceil(Tsw/Ts);
NTdown = ceil((Tsw-T)/(4*Ts));

t = 0:Ts:(M*NTsw*Ts-Ts);

m = B/(NT*Ts);
m_down = -B/(NTdown*Ts);

t_up = 0:Ts:((NT-1)*Ts);
t_down = 0:Ts:((NTdown-1)*Ts);

x = zeros(NTsw*M,1);
for l = 1:M
    x(((l-1)*NTsw+1):(NT+(l-1)*NTsw)) = sin(2*pi*t(((l-1)*NTsw+1):(NT+(l-1)*NTsw))+pi*m*(t_up.^2));
    x((NT+1+(l-1)*NTsw):(NT+NTdown+(l-1)*NTsw)) = sin(2*pi*t((NT+1+(l-1)*NTsw):(NT+NTdown+(l-1)*NTsw))+pi*m_down*(t_down.^2)+2*B*pi);
    x((NT+NTdown+1+(l-1)*NTsw):(NTsw*l)) = zeros(NTsw-NT-NTdown,1);
end

% x_linear = reshape(x,NTsw*M,1);
x_phase = angle(x);
x_freq = diff(x_phase)/(Ts*2*pi);
plot(t(1:(end)),x)