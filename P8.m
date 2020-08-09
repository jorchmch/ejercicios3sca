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
R=0.1*eye(2,2);
[G, ~, ~]=lqr(sys,Q,R)

Q1=eye(4,4);
R1=[0.01 0 ; 0 0.01];
[H, ~, ~] = lqe(sys,Q1,R1)


%cONTROLADOR k(S)

K11=[a-b*G];  % 4X4
K12=[-b*G];    %4X4
K21=[H*c]; %4x4
K22=[a-H*c-b*G];   %2x2

AK=[K11 K12;K21 K22];   %8x8
BK=[zeros(4,2);-H];     %8x2
CK=[c zeros(2,4)];     %2x8
DK=[0 0;0 0];
% % % % % % % % % % % % % % % 
sysK=ss(AK,BK,CK,DK)

