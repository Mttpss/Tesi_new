function [y,x,f] = echoSingleInterferenceFMCW(d,v,fc,fci,fs,B,Bi,NT,NTsw,Ti,Tswi,M,Mi,offset)
%ECHOSINGLETARGETFMCW Summary of this function goes here
%   Detailed explanation goes here
c = 3e8;
lambda_ci = c/fci;
Ts = 1/fs;
m = B/(NT*Ts);
NTi = floor(Ti/Ts);
NTswi = floor(Tswi/Ts);
mi = Bi/(NTi*Ts);

t = 0:Ts:(Mi*NTswi*Ts-Ts);
t_up = (0:Ts:((NTi-1)*Ts))';
R = d + v*t;
L = sqrt(10.^(-fspl(2*R,lambda_ci)/20)); % Free space path loss of the echo
t = reshape(t,NTswi,Mi);
L = reshape(L,NTswi,Mi);
R = reshape(R,NTswi,Mi);

x_tobeat = zeros(NTsw,M);
for l = 1:Mi
    x_tobeat(1:NT,l) = exp(1i*2*pi*(fc*t_up + m/2*(t_up.^2)))/sqrt(2);
    x_tobeat((NT+1):end,l) = zeros(NTsw-NT,1);
end
x_tobeat = reshape(x_tobeat,NTsw*M,1);

x = zeros(NTswi,Mi);
f = zeros(NTswi,Mi);
for l = 1:Mi
    x(1:NTi,l) = exp(1i*2*pi*(fci*t_up + mi/2*(t_up.^2)))/sqrt(2);
    x((NTi+1):end,l) = zeros(NTswi-NTi,1);
    f(1:NTi,l) = fci + mi/2*t_up;
    f((NTi+1):end,l) = zeros(NTswi-NTi,1);
end
x = reshape(x,NTswi*Mi,1);
disp(max(size(x)))
f = reshape(f,NTswi*Mi,1);
if offset < 0 
    x = [x(abs(offset)+1:end) ; zeros(abs(offset),1)];
    f = [f(abs(offset)+1:end) ; zeros(abs(offset),1)];
elseif offset > 0
    x = [zeros(abs(offset),1) ; x(1:(NTsw*M-abs(offset)))];
    f = [zeros(abs(offset),1) ; f(1:(NTsw*M-abs(offset)))];
end
disp(max(size(x)))
x = x(1:NTsw*M,1);
f = f(1:NTsw*M,1);

y_prebeat = zeros(NTswi,Mi);
for l = 1:Mi
    y_prebeat(1:NTi,l) = L(1:NTi,l).*exp(1i*2*pi*(- fci*d/c - fci*v*(l-1)*NTswi*Ts/c + (fci - mi*d/c - fci*v/c - mi*v*(l-1)*NTswi*Ts/c)*t_up + (mi/2 - mi*v/c)*(t_up.^2)))/sqrt(2);
    y_prebeat((NTi+1):end,l) = zeros(NTswi-NTi,1);
end
y_prebeat = reshape(y_prebeat,NTswi*Mi,1);
if offset < 0 
    y_prebeat = [y_prebeat(abs(offset)+1:end) ; zeros(abs(offset),1)];
elseif offset > 0
    y_prebeat = [zeros(abs(offset),1) ; y_prebeat(1:(NTsw*M-abs(offset)))];
end

y = x_tobeat.*conj(y_prebeat);
y = y(1:NTswi*Mi,1);
y = reshape(y,NTswi,Mi);
y = y(1:NTi,:);

end

