clc,clear;

% MIMO system
% define SISO systems
sys11=tf([-86.41],[15 0.37]);
sys12=tf([-0.05445],[1 0.524667 0.0873333 0.00185]);
sys21=tf([29324.85 327.75],[15 0.37]);
sys22=tf([0.25535],[1 0.524667 0.0873333 0.00185]);

% compose MIMO system
sys=[sys11 sys12;
     sys21 sys22];
[a,b,c,d]=tf2ss(sys);

%%%%%%%%% sensibilidad

al=[zeros(2,2) zeros(2,4); b a];
bl=[eye(2,2) ;zeros(4,2)];
cl=[zeros(2,2) c];
dl=eye(2,2);

w=logspace(-2,3,100);
%%% sensibilidad
sv = sigma(ss(al-bl*cl, bl, -cl,dl),w);
sv = 20*log10(sv);
figure
semilogx(w, sv)
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


sv = sigma(ss(al-bl*cl, bl, -cl, 0*dl),w);
sv = 20*log10(sv);
figure
semilogx(w, sv)

title('Complementary Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

