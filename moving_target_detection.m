close all;
clear all;
clc;

c=2.9979*10^8;  % Propagation velocity
fo=4.5e9;       % Carrier frequency (4.5GHz)
landa=c/fo;     % Wavelength (60cm for fo = 4.5e9)
theta_inc=pi/4; % Loook angle
Vplat=200;         % Velocity of platform
Rcenter=2000;  % Range distance to center of target area
l=0.5;           % Antenna vertical length actual
L=0.6;             %Antenna azimuth length

Swath=(landa*Rcenter)/(l*cos(theta_inc));

Rnear=Rcenter-Swath/2;     % Range Near
Rfar=Rcenter+Swath/2;   % Rang Far
Tr=(2*Swath)/(c*9);     %Pulse duration=>window length is 10 Tr
Ts=(2*Rnear)/c-Tr/2;             % Start time of sampling
Tf=(2*Rfar)/c+Tr/2;           % End time of sampling

Bw=100e6;%therefore resolution is 1.5m
Kr=Bw/Tr;%Modulation Rate

R0=Rcenter;
Timg=(landa*R0)/(L*Vplat);


max_PRF=(l*c)/(2*Rcenter*landa*tan(theta_inc));
min_PRF=2*Vplat/L;
%PRF=max_PRF+min_PRF/2;
PRF=1000;
Dta=1/PRF;
eta=-Timg/2:Dta:Timg/2;
Nffta=1024;
f=linspace(-PRF/2,PRF/2,Nffta);

Va=10; %along track velocity(m/s)
Vc=10; %cross track velocity(m/s)
Aa=0.1; %along track acceleration(m/s^2)
Ac=0.1; %cross track acceleration(m/s^2)
X=R0*sin(theta_inc);
Vr=X*Vc/R0;
Ar=X*Ac/R0;

tr=2*R0/c;

aux_var=((Vplat-Va)^2-Ar*R0) /(2*R0) ;
modulation_rate=2*((Vplat-Va)^2-Ar*R0) /(landa*R0) ;
dopler_centroid_shift=2*Vr/landa;
required_Bandwidth=dopler_centroid_shift+modulation_rate*Timg/2;
PRF_must_be=required_Bandwidth*2;

s1_aux=exp(1j*4*pi*Vr*eta/landa);
s2_aux=exp(-1j*4*pi*aux_var*eta.^2/landa);
s3_aux=exp(1j*4*pi*Kr*Vr*tr*eta/c);
s4_aux=exp(-1j*4*pi*Kr*aux_var*tr*eta.^2/c);
s5_aux=exp(-1j*4*pi*R0/landa);

%s1=s1_aux .* s2_aux .* s3_aux .* s4_aux .* s5_aux;
s1= s2_aux ;
plot(real(s1));
figure;
plot(f,abs(fftshift(fft(s1,Nffta))));


