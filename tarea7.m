clear 
clc 
close all

syms t s

num=[1 0 -2562 0 0];
den=[1 18 102.25 248 266.25 112.50];

G=tf(num,den);

To=0.01553;

[A B C D]=tf2ss(num,den);

%------ discretizacion truncada ------%
suma=0;
M=rank(A);
i=A*inv(A);

for n=1:1:M 
    s = (A^n*t^n)/n + suma;
    suma=s;
    Phi=i+s;
    Bi=Phi*B;
    gamma=int(Bi,0,t);
end

Ad=double(subs(Phi,t,To));
Bd=double(subs(gamma,t,To));

sys=ss(Ad,Bd,C,D,To);
Gz=tf(sys);
[numz,denz]=tfdata(Gz,'v');
%--------------------------------------------%
%----------- FCControlada -------------------%
Co=[B A*B A^2*B A^3*B A^4*B];
Gc=ss(Co,B,C,D);
GCo=tf(Gc);
% [numco,denco]=tfdata(GCo,'v');
%----------- Discretizacion -----------------%
M1=rank(Co);
i1=Co*inv(Co);

for n=1:1:M1 
    s1 = (Co^n*t^n)/n + suma;
    suma=s1;
    Phi1=i1+s1;
    Bi1=Phi1*B;
    gamma1=int(Bi1,0,t);
end

Ad1=double(subs(Phi1,t,To));
Bd1=double(subs(gamma1,t,To));

sys1=ss(Ad1,Bd1,C,D,To);
Gz1=tf(sys1);
[numco,denco]=tfdata(Gz1,'v');
%--------------------------------------------%
%----------- FCObservable -------------------%
Ob=[C;C*A;C*A^2;C*A^3;C*A^4];
Go=ss(Ob,B,C,D);
GOb=tf(Go);
%----------- Discretizacion -----------------%
M2=rank(Ob);
i2=Ob*inv(Ob);

for n=1:1:M2 
    s2 = (Ob^n*t^n)/n + suma;
    suma=s2;
    Phi2=i2+s2;
    Bi2=Phi2*B;
    gamma2=int(Bi2,0,t);
end

Ad2=double(subs(Phi2,t,To));
Bd2=double(subs(gamma2,t,To));

sys2=ss(Ad2,Bd2,C,D,To);
Gz2=tf(sys2);
[numob,denob]=tfdata(Gz2,'v');

figure(1)
y1=impz(numz,denz);
stem(y1)
title('Respuesta Impulso del sistema original');
grid on
figure(2)
y2=impz(numco,denco);
stem(y2)
title('Respuesta Impulso del sistema con FCControla');
grid on
figure(3)
y3=impz(numob,denob);
stem(y3)
title('Respuesta Impulso del sistema con FCObservable');
grid on