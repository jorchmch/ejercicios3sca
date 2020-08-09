// CONTROL AEREOGENERADOR
// Alumno: Jorch Mendoza Choquemamani
// Author: juan C. Cutipa-Luque
// Docente: juan C. Cutipa-Luque

//Codigo Reutilizado de un ejemplo Mostrado en Clase
clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows

A=    [0.0233    0.0052   -0.0650   -0.0280
   -0.1054    0.0014    0.1447    0.1004
         0   -0.3516   -0.0409   -0.3632
         0   -0.4921   -0.0228   -0.5331]
 
 
B=[    2.1926    0.3717
    1.9350    0.0000
    2.6300    0.0000
   -1.9909   -0.0000];

C=[         0         0   -2.8096   -0.8180
         0         0         0   13.2464];

D=[           0           0
        1955           0];
ap_s=A;  
bp_s=B;
cp_s=C;
dp_s=D;

// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap_s,bp_s)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap_s, cp_s)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

/*             Plot singular values of LTI the model                      */
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
poles=spec(ap_s)
tr  = trzeros(G)

w = logspace(-3,3);
sv = svplot(G,w);

scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


//g=minreal(g)
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
w = logspace(-3,3);
sv = svplot(G,w);

scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


ms=1.4;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4

//           --------     WS     ------------

s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
//Ws=syslin('c',ws)
Ws=blockdiag(ws1,ws2)

//           --------     WT     ------------

s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
Wt=syslin('c',wt)
//Wt=blockdiag(wt1,wt2)


//           --------     WR     ------------
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]



// ------------------ Plot weighting functions

svs = svplot(Ws,w);
svt = svplot(Wt,w);
scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");


[P,r]=augment(G,'ST'); // SENSIBILIDAD 
//[P,r]=augment(g,'SRT');
P = blockdiag(Ws,Wt,eye(G))*P;
//P=minreal(P);

//trick to tackle when "D12 is not full rank"
P.D(1,3)=0.0001;
P.D(2,4)=0.0001;

r=[2,2]
romin=0.0001
romax=2000
nmax=100;
//[K,ro]=h_inf(P,r,romin,romax,nmax)
// alternatives
//[AK,BK,CK,DK,(RCOND)] = hinf(P.A,P.B,P.C,P.D,2,2,4)
K = ccontrg(P, [2,2], 1.3) // this is good for me and for this system

// -------------- Analysis of the Feeedback Control System

[Se,Re,Te]=sensi(G,K) // S=(I+GK)^-1,    :sensibilidad
                      //T=I-S=GK(I+GK)^-1  :sensibilidad complementaria

// ------------------ Plot weighting functions

svS = svplot(Se,w);
svT = svplot(Te,w);
scf(3);
plot2d("ln", w,[-20*log(svs')/log(10) 20*log(svS')/log(10)],[-1 -1 2 2],leg="$\overline{\sigma}(W_S^{-1})$@$\underline{\sigma}(W_S^{-1})$@$\overline{\sigma}(S)$@$\underline{\sigma}(S)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
//set(gca(),"auto_clear","off")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");


scf(4);
plot2d("ln", w,[-20*log(svt')/log(10) 20*log(svT')/log(10)],[-1 -1 2 2],leg="$\overline{\sigma}(W_T^{-1})$@$\underline{\sigma}(W_T^{-1})$@$\overline{\sigma}(T)$@$\underline{\sigma}(T)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
xtitle("","Frequency (rad/s)", "Amplitude (dB)");

// --------------- Open loop and Closed loop analysis --------------

sysOL=G*K
sysCL=G*K*inv(eye(2,2)+G*K)
// eigenvalues in LHP (Left Half Plane), stable according to RH
spec(sysCL.A)

sv1= svplot(sysOL,w);
sv2= svplot(sysCL,w);
scf(5);
plot2d("ln", w,[20*log(sv1')/log(10) 20*log(sv2')/log(10)],[2 2 3 3],leg="$\overline{\sigma}(GK)$@$\underline{\sigma}(GK)$@$\overline{\sigma}(T)$@$\underline{\sigma}(T)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
//set(gca(),"auto_clear","off")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");

// Time responses in XCOS
