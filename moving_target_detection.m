%%
%%Doppler Domain Filtering
%%%written by A.Vafaei
close all;
clear all;
clc;
%%
%%%%%%Main Parameter
c=2.9979*10^8;  % Propagation velocity
fo=4.5e9;       % Carrier frequency (4.5GHz)
landa=c/fo;     % Wavelength (60cm for fo = 4.5e9)
theta_inc=pi/4; % Loook angle
Vplat=200;         % Velocity of platform
Rcenter=2000;  % Range distance to center of target area
l=0.5;           % Antenna vertical length actual
L=0.6;             %Antenna azimuth length,it is reduced to 0.6 because for estimation modulation rate via frft we need to apply bigger window therefor bigger Timg
%%
%%%%%%geometry of scene
Swath=(landa*Rcenter)/(l*cos(theta_inc));
Rnear=Rcenter-Swath/2;     % Range Near
Rfar=Rcenter+Swath/2;   % Rang Far
%%
%%%%%range time domain declaration
Tr=(2*Swath)/(c*9);     %Pulse duration=>window length is 10 Tr
Ts=(2*Rnear)/c-Tr/2;             % Start time of sampling
Tf=(2*Rfar)/c+Tr/2;           % End time of sampling
Bw=100e6;%therefore resolution is 1.5m
Kr=Bw/Tr;%Modulation Rate

R0=Rcenter;
tr=2*R0/c;

%%
%%%%%%azimuth time domain delaration
Timg=(landa*R0)/(L*Vplat);
max_PRF=(l*c)/(2*Rcenter*landa*tan(theta_inc));
min_PRF=2*Vplat/L;
%PRF=max_PRF+min_PRF/2;
PRF=2500;
Dta=1/PRF;
eta=-Timg/2:Dta:Timg/2;
Nffta=1024;
f=linspace(-PRF/2,PRF/2,Nffta);
%%
%%%%%%%%%%%%%%motion parameter for moving target

Va=10; %along track velocity(m/s)
Vc=10; %cross track velocity(m/s)
Aa=0.1; %along track acceleration(m/s^2)
Ac=0.1; %cross track acceleration(m/s^2)
X=R0*sin(theta_inc);
Vr=X*Vc/R0;
Ar=X*Ac/R0;


%%
%%%%%%%%%%%PRF 
aux_var=((Vplat-Va)^2-Ar*R0) /(2*R0) ;
modulation_rate=2*((Vplat-Va)^2-Ar*R0) /(landa*R0) ;
dopler_centroid_shift=2*Vr/landa;
required_Bandwidth1=dopler_centroid_shift+modulation_rate*Timg/2;
dopler_centroid_shift_rwm=2*Kr*Vr*tr/c;
modulation_rate_rcm=2*Kr*tr*((Vplat-Va)^2-Ar*R0) /(c*R0) ;
required_Bandwidth2=dopler_centroid_shift_rwm+modulation_rate_rcm*Timg/2;
required_Bandwidth=required_Bandwidth1+required_Bandwidth2;
%required_Bandwidth=dopler_centroid_shift+dopler_centroid_shift_rwm+max(modulation_rate,modulation_rate_rcm)*Timg/2;
PRF_must_be=required_Bandwidth*2;


%%
%%%%%%%%signal of moving target
s1_aux=exp(1j*4*pi*Vr*eta/landa);
s2_aux=exp(-1j*4*pi*aux_var*eta.^2/landa);
s3_aux=exp(1j*4*pi*Kr*Vr*tr*eta/c);
s4_aux=exp(-1j*4*pi*Kr*aux_var*tr*eta.^2/c);
s5_aux=exp(-1j*4*pi*R0/landa);

s1=s1_aux .* s2_aux .* s3_aux .* s4_aux .* s5_aux;
%s1= s2_aux;
figure;
plot(real(s1));
figure;
plot(f,abs(fftshift(fft(s1,Nffta))));
%%
%motion parameters are set to zero so signal model for moving Target can
%be used for stationary target

Va=0; %along track velocity(m/s)
Vc=0; %cross track velocity(m/s)
Aa=0; %along track acceleration(m/s^2)
Ac=0; %cross track acceleration(m/s^2)
X=R0*sin(theta_inc);
Vr=X*Vc/R0;
Ar=X*Ac/R0;
aux_var=((Vplat-Va)^2-Ar*R0) /(2*R0) ;

s1_aux=exp(1j*4*pi*Vr*eta/landa);
s2_aux=exp(-1j*4*pi*aux_var*eta.^2/landa);
s3_aux=exp(1j*4*pi*Kr*Vr*tr*eta/c);
s4_aux=exp(-1j*4*pi*Kr*aux_var*tr*eta.^2/c);
s5_aux=exp(-1j*4*pi*R0/landa);

%s2=s1_aux .* s2_aux .* s3_aux .* s4_aux .* s5_aux;
s2= s2_aux;
figure;
plot(real(s2));
figure;
plot(f,abs(fftshift(fft(s2,Nffta))));

%%
%main signal

s=s1+s2;
%s1= s2_aux
figure;
plot(real(s));
figure;
plot(f,abs(fftshift(fft(s,Nffta))));

%%
%%%%%%%filter design

fc=400;
fs=PRF;
omega=2*pi*fc/fs;
index=round(fc*Nffta/fs);
H_f_temp=ones(1,Nffta/2);
H_f_temp(1:index)=zeros(1,index);
H_f_temp2=(flipud(H_f_temp'))';
H_f_filter=[H_f_temp2,H_f_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%apply filter to signals
s1_f=fftshift(fft(s1,Nffta));
s2_f=fftshift(fft(s2,Nffta));
s_f=fftshift(fft(s,Nffta));
s1_f_after_filter=s1_f .* H_f_filter;
s2_f_after_filter=s2_f .* H_f_filter;
s_f_after_filter=s_f .* H_f_filter;
figure;
plot(f,abs(s1_f_after_filter));
ylim([0 80]);
figure;
plot(f,abs(s2_f_after_filter));
ylim([0 80]);
figure;
plot(f,abs(s_f_after_filter));
ylim([0 80]);

%%
%%%%%%%%%%%%fractional fourier Transform

modulation_rate_stationary_target=2*((Vplat-Va)^2-Ar*R0) /(landa*R0) ;
a=2*(atan(modulation_rate_stationary_target/fs)+pi/2)/pi;
s2_after_frft=frft(s2,a);
plot(abs(s2_after_frft));

%%
%%%%%%%%%%%%%%%%test
alpha=-2:0.001:2;
max_temp=0;
angle=0;
counter=1;

for a=alpha
    Y=frft(s2,a);
    if(max_temp<max(abs(Y)))
       angle=a;
       max_temp=max(abs(Y));
    end
end
Y=frft(s2,angle);
figure;
plot(abs(Y));
figure;
plot(abs(fftshift(fft(s2))));
phi=angle*pi/2;
kr_estimate=tan(phi-pi/2)*fs;








