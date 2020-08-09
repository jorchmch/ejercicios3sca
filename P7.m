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
%Q =c'*c;
Q=eye(4,4);
R=[0.01 0 ; 0 0.01];
[H, ~, ~] = lqe(sys,Q,R)

