clc,clear;


w=logspace(-1,3,400);
a=70;b=130;c=200;
x1=20*ones(1,a); x2=60*zeros(1,b); x3=-20*ones(1,c);
xt=[x1 x2 x3];
x4=[5*ones(1,100) 0*zeros(1,300)];
% barreras de estabilidad

semilogx(w, xt,'r')
ylim([-60,60])
grid;
hold on;
semilogx(w,x4,'b')
title('Barreras de Estabilidad')
xlabel('w(rad/s)')
ylabel('dB')

hold off;
