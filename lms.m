clc
close all
clear all

N= 296; %input('length of sequence N = ');
t=[0:N-1];
w0=0.001;  phi=0.1;
d=sin(2*pi*[1:N]*w0+phi);
x=d+randn(1,N)*0.5;
w=zeros(1,N); 
mu= 0.06 ; %0.0005; %input('mu = ');
for i=1:N
   e(i) = d(i) - w(i)' * x(i);
   w(i+1) = w(i) + mu * e(i) * x(i);
end
for i=1:N
yd(i) = sum(w(i)' * x(i));  
end

d_est = w(1,N) * x ;

subplot(221),plot(t,d),ylabel('Desired Signal'),
subplot(222),plot(t,x),ylabel('Input Signal+Noise'),
subplot(223),plot(t,e),ylabel('Error'),
subplot(224),plot(t,yd),ylabel('Adaptive Desired output');

figure();
plot(t,d_est),ylabel('Estimated output');

