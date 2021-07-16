%BP_Chebyshev

clear
clc

%Τρέχοντας το σκριπτάκι παίρνουμε όλα τα επιθυμητά πλοταρίσματα
%-Τις ΣΜ κάθε μονάδας καθώς και της ολικής
%-Όλες τις επιμέρους ΣΜ σε κοινό πλοτ
%-Συνάρτηση απόσβεσης της ολικής ΣΜ
%-Σήμα εισόδου και εξόδου σε κοινό πλοτ
%-Fourier σήματος εισόδου
%-Fourier σήματος εξόδου

%Προδιαγραφές
a1 = 8;
a2 = 9;
a3 = 8;
a4 = 1;

f0 = 900;
f1 = 650+25*a3;
f2 = f0^2/f1;
D = 2.2*(f0^2 -f1^2)/f1;
f3 = (-D+sqrt(D^2+4*f0^2))/2;
f4 = f0^2/f3;

%Κυκλικές συχνότητες
w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

amin = 28.5+a4*5/9;
amax = 0.5+a3/36;


W_p = 1;
W_s = (w4 - w3)/(w2-w1);

%Υπολογισμός τάξης φίλτρου
n = acosh(((10^(amin/10)-1)/(10^(amax/10)-1))^(1/2))/acosh(W_s);
n = ceil(n);

e = sqrt(10^(amax/10)-1);

%Συχνότητας ημίσειας ισχύος και συντελεστή α
W_hp = cosh(acosh(1/e)/n);
a = asinh(1/e)/n;

%Γωνίες Βutterworth 4ης τάξης
ps_k1 = 22.5;
ps_k2 = -22.5;
ps_k3 = 67.5;
ps_k4 = -67.5;

%Πόλοι Chebyshev
s_k1 = -sinh(a)*cosd(ps_k1) + j*cosh(a)*sind(ps_k1);
s_k2 = -sinh(a)*cosd(ps_k2) + j*cosh(a)*sind(ps_k2);
s_k3 = -sinh(a)*cosd(ps_k3) + j*cosh(a)*sind(ps_k3);
s_k4 = -sinh(a)*cosd(ps_k4) + j*cosh(a)*sind(ps_k4);

%Μετασχηματισμός πόλων
bw = w2-w1;
qc = w0/bw;

%Μετασχηματισμός του πόλου s1,2
S2 = abs(real(s_k1));
OMEGA2 = imag(s_k1);
C2 = (S2)^2+(OMEGA2)^2;
D2 = 2*S2/qc;
E2 = 4+C2/((qc)^2);
G2 = sqrt(E2^2 -4*D2^2);
Q_S2 = sqrt((E2+G2)/2)/D2;
k2 = S2*Q_S2/qc;
W2 = k2+sqrt(k2^2-1);
w_02 = W2*w0;
w_01 =w0/W2;

%Μετασχηματισμός του πόλου s3,4
S4 = abs(real(s_k4));
OMEGA4 = imag(s_k4);
C4 = (S4)^2+(OMEGA4)^2;
D4 = 2*S4/qc;
E4 = 4+C4/((qc)^2);
G4 = sqrt(E4^2 -4*D4^2);
Q_S4 = sqrt((E4+G4)/2)/D4;
k4 = S4*Q_S4/qc;
W4 = k4+sqrt(k4^2-1);
w_04 = W4*w0;
w_03 =w0/W4;

%Aφού Q_S2>5 και Q_S4>5 χρησιμοποιούμε την τεχνική Q-enhancement για όλες
%τις Μονάδες

%1η Μονάδα                
b1 = 1;                   %    0.01<b<100      
R11_old = 1/(sqrt(b1));
R12_old = sqrt(b1);
k1 = (Q_S2*(b1+2)-sqrt(b1))/(2*Q_S2-sqrt(b1));
C11_old = 1;
C12_old = 1;
RA1_old = 1;
RB1_old = (k1-1)*RA1_old;
% Κλιμακοποίηση
kf1  = w_01;
km1  = C11_old/(kf1*0.01*10^(-6));
C11_new = C11_old/(kf1*km1);
C12_new = C12_old/(kf1*km1);
R11_new = km1 * R11_old;
R12_new = km1 * R12_old;
RA1_new  = km1*RA1_old;
RB1_new  = km1*RB1_old;


%2η Μονάδα                      
b2   = 2;                      
R21_old = 1/(sqrt(b2));
R22_old = sqrt(b2);
k2   = (Q_S2*(b2+2)-sqrt(b2))/(2*Q_S2-sqrt(b2));
C21_old = 1;
C22_old = 1;
RA2_old = 1;
RB2_old = (k2-1)*RA2_old;
% Κλιμακοποίηση
kf2  = w_02;
km2  = C21_old/(kf2*0.01*10^(-6));
C21_new = C21_old/(kf2*km2);
C22_new = C22_old/(kf2*km2);
R21_new = km2 * R21_old;
R22_new = km2 * R22_old;
RA2_new  = km2*RA2_old;
RB2_new  = km2*RB2_old;



%3η Μονάδα                   
b3  = 2;                   
R31_old = 1/(sqrt(b3));
R32_old = sqrt(b3);
k3   = (Q_S4*(b3+2)-sqrt(b3))/(2*Q_S4-sqrt(b3));
C31_old = 1;
C32_old = 1;
RA3_old = 1;
RB3_old = (k3-1)*RA3_old;

%Κλιμακοποίηση
kf3  = w_03;
km3  = C31_old/(kf3*0.01*10^(-6));
C31_new = C31_old/(kf3*km3);
C32_new = C32_old/(kf3*km3);
R31_new = km3 * R31_old;
R32_new = km3 * R32_old;
RA3_new  = km3*RA3_old;
RB3_new  = km3*RB3_old;


%4η Μονάδα                      
b4 = 2;                 
R41_old = 1/(sqrt(b4));
R42_old = sqrt(b4);
k4   = (Q_S4*(b4+2)-sqrt(b4))/(2*Q_S4-sqrt(b4));
C41_old = 1;
C42_old = 1;
RA4_old = 1;
RB4_old = (k4-1)*RA4_old;

%Κλιμακοποίηση
kf4  = w_04;
km4  = C41_old/(kf4*0.01*10^(-6));
C41_new = C41_old/(kf4*km4);
C42_new = C42_old/(kf4*km4);
R41_new = km4 * R41_old;
R42_new = km4 * R42_old;
RA4_new  = km4*RA4_old;
RB4_new  = km4*RB4_old;


%Κέρδοι
H1   = k1*b1/(2*(k1-1)-b1);
H2   = k2*b2/(2*(k2-1)-b2);
H3   = k3*b3/(2*(k3-1)-b3);
H4   = k4*b4/(2*(k4-1)-b4);


%Συναρτήσεις Μεταφοράς Moνάδων 
TF_1 = tf ([0 H1*(w_01 /Q_S2) 0], [1, w_01/Q_S2, w_01^2]);
TF_2 = tf ([0 H2*(w_02 /Q_S2) 0], [1, w_02/Q_S2, w_02^2]);
TF_3 = tf ([0 H3*(w_03 /Q_S4) 0], [1, w_03/Q_S4, w_03^2]);
TF_4 = tf ([0 H4*(w_04 /Q_S4) 0], [1, w_04/Q_S4, w_04^2]);



%Ρύθμιση Κέρδους
h_TOT = bode(TF_1,w0)*bode(TF_2,w0)*bode(TF_3,w0)*bode(TF_4,w0);	
a = 1 /h_TOT;
Z_a = R11_new/a;
Z_b = R11_new/(1-a);

%Ολική ΣΜ
T_Total =  a*TF_1 * TF_2 * TF_3*  TF_4;
%T_Total_before = TF_1 * TF_2 * TF_3*  TF_4

plot_transfer_function(TF_1, [ f0 f1 f2 f3 f4]);
plot_transfer_function(TF_2, [ f0 f1 f2 f3 f4]);
plot_transfer_function(TF_3, [ f0 f1 f2 f3 f4]);
plot_transfer_function(TF_4, [ f0 f1 f2 f3 f4]);
plot_transfer_function(T_Total, [ f0 f1 f2 f3 f4]);
plot_transfer_function(inv(T_Total), [ f0 f1 f2 f3 f4]);
%plot_transfer_function(inv(T_Total_before), [ f0 f1 f2 f3 f4]);΄
ltiview({'bodemag'}, TF_1, TF_2, TF_3, TF_4, T_Total);

%Fourier
Fs = 100000;
t = 0:1/Fs:0.03;
T = 1/2000;
%Δημιουργία σήματος εισόδου
u = cos((w0 - (w0-w1)*0.5)*t) + 0.8*cos((w0+(w0+w1)/3)*t) + 0.8*cos(0.4*w3*t)+0.6*cos(2.5*w4*t)+0.5*cos(3*w4*t);
N = length(u);
Fbins = ((0: 1/N: 1-1/N)*Fs);
U=abs(fft(u));
Vout = lsim(T_Total, u, t);
%Σήμα εισόδου και σήμα εξόδου
figure
plot(t,u)
hold on
plot(t,Vout)
xlabel('Time(s)');
ylabel('Voltage(V)');
hold off
%Fourier Εισόδου
figure
plot(Fbins,U)
xlabel('Frequency(Hz)'); 
ylabel('Magnitude');
F=abs(fft(Vout));
%Fourier Εξόδου
figure
plot(Fbins,F)
xlabel('Frequency(Hz)'); 
ylabel('Magnitude');












