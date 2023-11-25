x_forward = linspace(0,1,100);
x_backward = linspace(-3,0,100);
x_backward(end) = [];
x = [x_backward x_forward];
seno = 2*pi*5*x;

figure
plot(x,seno);