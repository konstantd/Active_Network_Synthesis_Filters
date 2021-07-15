%LP_Chebyshev

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
m=3;
fp = (3+m)* 1000;
fs = 1.72*fp;
amin = 24 + (max(1, a3)-5)*(3/4);
amax = 0.65 +(max(1, a4)-5)/10; 

e= sqrt(10^(amax / 10) - 1);

%Κυκλικές συχνότητες
wp= 2*pi*fp;
ws= 2*pi*fs;
WS= ws/wp;

%Υπολογισμός τάξης φίλτρου
n  = acosh( sqrt( ( 10^(amin / 10) - 1 ) / ( 10^(amax / 10) - 1) ) ) / acosh(WS);
n = ceil(n);

%Συχνότητας ημίσειας ισχύος και συντελεστή α
Whp = cosh( (1 / n) * acosh( 1 / e) );
a= asinh(1/e)/n;

%Η πραγματική ωhp
whp = 2 * pi * Whp * fp; 

%οι πόλοι Chebyshev
p1 = -sinh(a);
p2 = -sinh(a)*cosd(36) + j*cosh(a)*sind(36);
p3 = -sinh(a)*cosd(36) - j*cosh(a)*sind(36);
p4 = -sinh(a)*cosd(72) + j*cosh(a)*sind(72);
p5 = -sinh(a)*cosd(72) - j*cosh(a)*sind(72);

%Υπολογισμός Q της κάθε μονάδας
Q1  = abs(p1)/(2*abs(real(p1))); 
Q23 = abs(p2)/(2*abs(real(p2))); 
Q45 = abs(p4)/(2*abs(real(p4))); 

%Υπολογισμός Ω της κάθε μονάδας
OMEGA_1  = abs(p1);
OMEGA_23 = abs(p2);
OMEGA_45 = abs(p4);

%Πραγματικά μέτρα πόλων
OMEGA_01 = OMEGA_1 * wp;
OMEGA_02 = OMEGA_23 * wp;
OMEGA_03 = OMEGA_45 * wp;

%1η Μονάδα και κλιμακοποίηση της 
R11_old=1;
C11_old=1/(R11_old*abs(p1));
kf=wp;                             
km=C11_old/(kf*1e-06);
R11_new=km*R11_old;
C11_new=C11_old/(km*kf);


%2η Μονάδα  και κλιμακοποίηση της 
C21_old= 2*Q23;
C22_old= 1/(2*Q23);
R21_old= 1;
R22_old= 1;
kf = wp*OMEGA_23;
km= C21_old/(kf*1e-06);
C21_new=C21_old/(km*kf);
C22_new=C22_old/(km*kf);
R21_new=km;
R22_new=km;


%3η Μονάδα  και κλιμακοποίηση της 
C31_old=2*Q45;
C32_old=1/(2*Q45);
R31_old=1;
R32_old=1;
kf = wp*OMEGA_45;
km= C31_old/(kf*1e-06);
C31_new=C31_old/(km*kf);
C32_new=C32_old/(km*kf);
R31_new=km;
R32_new=km;

%Συναρτήσεις Μεταφοράς Μονάδων
T1 = tf( OMEGA_01, [ 1 OMEGA_01 ] );
T2 = tf( OMEGA_02 ^2 , [ 1 ( OMEGA_02 / Q23 ) OMEGA_02 ^2 ] );
T3 = tf( OMEGA_03 ^2 , [ 1 ( OMEGA_03 / Q45 ) OMEGA_03 ^2 ] );

%Oλική ΣΜ
TLP = series( T1, T2);
TLP = series( TLP, T3);

plot_transfer_function( T1, [fs whp/(2*pi) fp]  );
plot_transfer_function( T2, [fs whp/(2*pi) fp]  );
plot_transfer_function( T3, [fs whp/(2*pi) fp]  );
plot_transfer_function( TLP, [fs whp/(2*pi) fp] );
plot_transfer_function( inv(TLP), [fs whp/(2*pi) fp] );
ltiview({'bodemag'}, T1, T2, T3, TLP)


%Fourier
Fs = 100000;
t = 0:1/Fs:0.004;
T = 1/2000;
%Δημιουργία τετραγωνικού παλμού
pulsewidth = 0.4*T;
pulseperiods = [0:20]*T;
u = pulstran(t,pulseperiods,@rectpuls,pulsewidth);
N = length(u);
Fbins = ((0: 1/N: 1-1/N)*Fs);
U=abs(fft(u));
Vout = lsim(TLP, u, t);
%Σήμα εισόδου και σήμα εξόδου
figure
plot(t,u)
axis([0 0.004 -0.5 1.5])  
hold on
plot(t,Vout)
xlabel('Time(s)')
ylabel('Voltage(V)')
hold off
%Fourier Εισόδου
figure
plot(Fbins,U)
xlabel('Frequency(Hz)')
ylabel('Magnitude')
F=abs(fft(Vout));
%Fourier Εξόδου
figure
plot(Fbins,F)
xlabel('Frequency(Hz)')
ylabel('Magnitude')



