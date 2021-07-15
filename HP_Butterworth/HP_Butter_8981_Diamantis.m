%HP_Butterworth

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
m = 1;

fp = ( 4 + m ) * 1000;
fs = fp / 2.6;
amin = 24+ ( a3 * ( 6 / 9 ) );
amax = 0.5 + ( a4 / 36 );
ws = 2 * pi * fs;
wp = 2 * pi * fp;
Wp = 1;
Ws = wp / ws;

%Υπολογισμός τάξης φίλτρου
n  = log10 ( ( 10^( amin / 10 ) - 1 ) / ( 10^( amax / 10 ) - 1 ) ) / ( 2 * log10 (Ws) );
n = ceil(n);

%Συχνότητας ημίσειας ισχύος 
Whp = 1 / ( ( 10^ ( amax / 10 ) - 1 )^( 1 / ( 2 * n ) ) );
wo = wp/ Whp;


%Πόλοι
p1= -1;
p2 = -0.809 + 0.587 * 1i;
p3 = -0.809 - 0.587 * 1i;
p4 = -0.309 + 0.951 * 1i;
p5 = -0.309 - 0.951 * 1i;

Q1=0.5;
Q23=0.628;
Q45= 1.628;


W1 = 1;
W23=1;
W45=1;


%1η Μονάδα και κλιμακοποίηση της 
C11_old=1;
R11_old=1;
kf=wo;
km= C11_old/( kf * 0.01* 10^(-6));
k1=1;
C11_new=1/(kf*km);
R11_new=km;

%2η Μονάδα και κλιμακοποίηση της 
C21_old = 1;
C22_old = 1;
R21_old = 1 / ( 2 * Q23);
R22_old = 2 * Q23;
kf2 = wo;
km2 =  km;
R21_new = km2 * R21_old;
R22_new = km2 * R22_old;
C21_new = C21_old / ( kf2 * km );
C22_new = C22_old / ( kf2 * km );

%3η Μονάδα και κλιμακοποίηση της 
C31_old = 1;
C32_old = 1;
R31_old= 1/(2 * Q45);
R32_old= 2 *  Q45;
km3= km;
R31_new= R31_old * km;
R32_new=R32_old * km;
C31_new = C31_old / ( kf2 * km );
C32_new = C32_old / ( kf2 * km );

%Ρύθμιση Κέρδους
k=1;
aGain = 10^( 10 / 20 ) / k;
Za = 1000;
Zb = Za * (aGain-1);


%Συναρτήσεις Μεταφοράς Moνάδων 
T1 = tf([1 0],[1 (wo)]);
T2 = tf([1 0 0],[1 ( wo / Q23 ) (wo)^2]);
T3 = tf([1 0 0],[1 ( wo / Q45 ) (wo)^2]);

%Oλική ΣΜ
TLP = aGain*T1*T2*T3;

plot_transfer_function( T1, [fs Whp/(2*pi) fp]  );
plot_transfer_function( T2, [fs wo/(2*pi) fp]  );
plot_transfer_function( T3, [fs wo/(2*pi) fp]  );
plot_transfer_function( TLP, [fs wo/(2*pi) fp] );
plot_transfer_function( inv(TLP), [fs wo/(2*pi) fp] );
ltiview({'bodemag'}, T1, T2, T3, TLP);



%Fourier
Fs = 100000;
t = 0:1/Fs:0.004;
T = 1/2000;
%Δημιουργία σήματος εισόδου
u=cos(0.2*ws*t)+0.6*cos(0.7*ws*t)+1.5*cos(1.6*wp*t)+0.7*cos(2.4*wp*t)+0.4*cos(3.5*wp*t);
N = length(u);
Fbins = ((0: 1/N: 1-1/N)*Fs);
U=abs(fft(u));
Vout = lsim(TLP, u, t);
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
xlabel('Frequency(Hz)') 
ylabel('Magnitude')




