%BE_Inverse_Chebyshev

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
a1=8;
a2=9;
a3=8;
a4=1;

f0=2.4*1000;
f1=1725 + 25*a4;
D= (1/2.2)*(f0^2-f1^2)/f1;
f2= f0^2/f1;
f3= (-D+ sqrt(D^2+4*f0^2))/2;
f4=f0^2/f3;
amin= 28+a3*5/9;
amax= 0.5+ a4/18;

%Κυκλικές συχνότητες
w1 = 2*pi*f1;
w2= 2*pi*f2;
w3= 2*pi*f3;
w4 = 2*pi*f4;
wo = 2*pi*f0;
Wp = 1;

Ws = ( w2 - w1 ) / ( w4 - w3 );
BW = 2*pi*(f2-f1);
qc = wo / BW;

%Υπολογισμός τάξης φίλτρου
n = acosh( sqrt( ( 10^( amin / 10 ) -1 ) / ( 10^( amax / 10 ) - 1 ) ) ) / acosh( Ws );
n = ceil(n);

%Συχνότητας ημίσειας ισχύος και συντελεστή α και ε
e = 1 / sqrt( 10^( amin / 10 ) -1 );
a = ( 1 / n ) * ( asinh( 1 / e ) );
whp = 1 / cosh( acosh( 1 / e ) / n);

%Γωνίες Βutterworth 4ης τάξης
ps_k1 = 22.5;
ps_k2 = -22.5;
ps_k3 = 67.5;
ps_k4 = -67.5;

%Πόλοι Chebyshev 
p1 = -sinh(a)*cosd(ps_k1) + j*cosh(a)*sind(ps_k1);
p2 = -sinh(a)*cosd(ps_k2) + j*cosh(a)*sind(ps_k2);
p3 = -sinh(a)*cosd(ps_k3) + j*cosh(a)*sind(ps_k3);
p4 = -sinh(a)*cosd(ps_k4) + j*cosh(a)*sind(ps_k4);

%Yπολογισμος Ω και Q
OMEGA_12 = abs(p1);
OMEGA_34 = abs(p3);

Q12 = OMEGA_12 /(2*abs(real(p1))); 
Q34 = OMEGA_34 / (2*abs(real(p3))); 

%Πόλοι ICΗ
InvW1 = 1 / OMEGA_12;
InvW3 = 1 / OMEGA_34;


%Κλιμακοποίηση πόλων 
InvW1= InvW1 * Ws;
InvW3= InvW3 * Ws;

%Μηδενικά
Z1 = sec( pi /(2*n));
Z2 = sec( 3*pi / (2*n));

%Κλιμακοποίηση μηδενικών
Z1 = Z1 * Ws;
Z2 = Z2 * Ws;

%Αντιστροφή πόλων 
InvW1 =  1 / (  InvW1 );
InvW3 =  1 / (  InvW3 );


%Θέσεις πόλων
S12= -InvW1 / (2* Q12);
W12 = sqrt( InvW1^2 - S12^2 );

S34= -InvW3 / (2* Q34);
W34 = sqrt( InvW3^2 - S34^2 );


%Αντιστροφή μηδενικών 
Z1= 1 / Z1;
Z2= 1 / Z2;


%Μετασχηματισμός του πόλου s1,2
p12Inv = S12 + ( W12 * 1i );
C1 = S12^2 + W12^2;
D1= -2* S12 / qc;
E1= 4 + C1/ qc^2;
G1= sqrt (E1^2 - 4* D1^2);
Q1_2= 1/D1 * sqrt ( 1/2* ( E1+ G1) );
k1= -S12 * Q1_2 /qc;
W1= k1 + sqrt( k1^2 -1);
wo1 = 1/W1 * wo;
wo2 = W1 * wo;

%Μετασχηματισμός του πόλου s3,4
p34Inv = S34 + ( W34 * 1i );
C2= S34^2 + W34^2;
D2= -2* S34 / qc;
E2= 4 + C2/(qc^2);
G2= sqrt (E2^2 - 4* D2^2);
Q3_4= 1/D2 * sqrt ( 1/2* ( E2+ G2) );
k2= -S34 * Q3_4 /qc;
W3_4= k2 + sqrt(k2^2 -1);
wo3 = 1/W3_4 * wo;
wo4 = W3_4 * wo;


%Μετασχηματισμός μηδενικών
%Z1
Kz1   = 2 + (Z1^2) / (qc^2);
x1    = ( Kz1 + sqrt( Kz1^2 - 4 ) ) / 2;
wz1 = wo * ( sqrt(x1) );
wz2 = wo / ( sqrt(x1) );

%Z2
Kz2   = 2 + (Z2^2) / (qc^2);
x2    = ( Kz2 + sqrt( Kz2^2 - 4 ) ) / 2;
wz3 = wo * ( sqrt(x2) );
wz4 = wo / (sqrt(x2));

%Έλεγχος για LPN ή HPN
wl1=wz1/wo1;
wl2=wz2/wo2;
wl3=wz3/wo3;
wl4=wz4/wo4;


%Moνάδα 1               wl1>1  -----> LPN 
R11_old=1;
R14_old=1;
R12_old= 4*Q1_2^2;
R13_old=  wl1^2 / (2* Q1_2^2);
R15_old=4*Q1_2^2/(wl1^2-1);
k11 =1/(R13_old +1); 
C11_old=1/(2*Q1_2);
H1 = k11*wl1^2;

%Κλιμακοποίηση
kf1=wo1;
km1=C11_old/(kf1*0.1e-06);
R11_new=km1*R11_old;
R12_new=R12_old*km1;
R13_new=km1*R13_old;
R14_new=R14_old*km1;
R15_new = R15_old*km1;
C11_new=C11_old/(km1*kf1);
C12_new=C11_new;

%Moνάδα 2             wl1<1  -----> HPN
k21 = 1/wl2^2 - 1;
k22 = ((2 + k21)*Q1_2^2) /((2 + k21)*Q1_2^2 + 1);
k23= k22/wl2^2;
R21_old = 1;
R22_old = ((2 + k21)^2)*Q1_2^2;
R23_old = 1;
R24_old = (2 + k21)*Q1_2^2;
C22_old = 1/((2+k21)*Q1_2);
C23_old = C22_old;
C21_old = k21*C22_old;
H2 = k23*wl2^2; 

%Κλιμακοποίηση
kf2 = wo2;
km2 = C22_old/(kf2* 0.1e-06);
R21_new = R21_old*km2;
R22_new= R22_old*km2;
R23_new = R23_old*km2;
R24_new = R24_old*km2;
C22_new = C22_old/(km2*kf2);
C23_new = C23_old/(km2*kf2);
C21_new = C21_old/(km2*kf2);


%Moνάδα 3               wl3>1  -----> LPN 
R31_old=1;
R34_old=1;
R32_old= 4*Q3_4^2;
R33_old=  wl3^2 / (2* Q3_4^2); 
R35_old=4*Q3_4^2/(wl3^2-1);
k31 =1/(R33_old +1);    
C31_old=1/(2*Q3_4);
H3 = k31*wl3^2;

%Κλιμακοποίηση
kf3=wo3;
km3=C31_old/(kf3*0.1e-06);
R31_new=km3*R31_old;
R32_new=R32_old*km3;
R33_new=km3*R33_old;
R34_new=R34_old*km3;
R35_new = R35_old*km3;
C31_new=C31_old/(km3*kf3);
C32_new=C31_new;


%Moνάδα 4           wl1<1  -----> HPN
k41 = 1/wl4^2 - 1;
k42 = (2 + k41)*Q3_4^2/((2 + k41)*Q3_4^2 + 1);
k43= k42/wl4^2; 
R41_old = 1;
R42_old = ((2 + k41)^2)*Q3_4^2;
R43_old = 1;
R44_old = (2 + k41)*Q3_4^2;
C42_old = 1/((2+k41)*Q3_4);
C43_old = C42_old;
C41_old = k41*C42_old;
H4 = k43  *wl4^2; 

%Κλιμακοποίηση
kf4 = wo4;
km4 = C42_old/(kf4* 0.1e-06);
R41_new = R41_old*km4;
R42_new= R42_old*km4;
R43_new = R43_old*km4;
R44_new = R44_old*km4;
C42_new = C42_old/(km4*kf4);
C43_new = C43_old/(km4*kf4);
C41_new = C41_old/(km4*kf4);


%Συναρτήσεις Μεταφοράς Moνάδων 
T1 = tf( [H1 0 ( H1 * wz1^2 ) ], [ 1 ( wo1 / Q1_2 ) wo1^2 ] );
T2 = tf( [H2 0 ( H2 * wz2^2 ) ], [ 1 ( wo2 / Q1_2 ) wo2^2 ] );
T3 = tf( [H3 0 ( H3 * wz3^2 ) ], [ 1 ( wo3 / Q3_4 ) wo3^2 ] );
T4 = tf( [H4 0 ( H4 * wz4^2 ) ], [ 1 ( wo4 / Q3_4 ) wo4^2 ] );  

TBE = T1*T2*T3*T4;

%Ρύθμιση Κέρδους
K = H1*H2*H3*H4;
again = (10^(0.5))/K;   
Rb = 100;
Ra = again*Rb-Rb;

%Ολική ΣΜ
TBE = again * TBE;


plot_transfer_function( T1, [f1 f3 f0 f4 f2] );
plot_transfer_function( T2, [f1 f3 f0 f4 f2] );
plot_transfer_function( T3, [f1 f3 f0 f4 f2] );
plot_transfer_function( T4, [f1 f3 f0 f4 f2] );
plot_transfer_function( TBE, [f1 f3 f0 f4 f2] );
ltiview({'bodemag'}, T1, T2, T3,T4, TBE)
plot_transfer_function( inv(TBE), [f1 f3 f0 f4 f2] );


%Fourier
Fs = 100000;
t = 0:1/Fs:0.004;
T = 1/2000;
%Δημιουργία σήματος εισόδου
u= 0.8*cos((wo-(wo-w3)/2)*t)+cos((wo+(wo+w3)/2)*t)+ cos(0.5*w1*t)+ 0.8*cos(2.4*w2*t)+ 0.4*cos(3.5*w2*t);
N = length(u);
Fbins = ((0: 1/N: 1-1/N)*Fs);
U=abs(fft(u));
Vout = lsim(TBE, u, t);
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




