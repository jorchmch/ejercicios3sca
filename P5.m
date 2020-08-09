% MIMO system responses

% MIMO system
% define SISO systems
sys11=tf([-86.41],[15 0.37]);
sys12=tf([-0.05445],[1 0.524667 0.0873333 0.00185]);
sys21=tf([29324.85 327.75],[15 0.37]);
sys22=tf([0.25535],[1 0.524667 0.0873333 0.00185]);

% compose MIMO system
sys=[sys11 sys12;
     sys21 sys22];

% compute singular values
% and plot
w=logspace(-2,2,100);
sv=sigma(sys,w);
V=sv(:,51) %VALORES SINGULARES EN 1 RAD/S


sv=20*log10(sv);
figure;
semilogx(w,sv);
grid

