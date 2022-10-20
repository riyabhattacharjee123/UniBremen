%% simple LMS %%
%simple_lms.m%

close all;
clear all;
clc;

t = 0.001:0.001:1;
D = 2*sin(2*pi*50*t);

n = size(D,2);
A = D(1:n)+ 0.9*randn(1,n);

M = 25;

w = zeros(1,M);
wi = zeros(1,M);

E=[];
mu = 0.0005;

for i = M:n
    j =  A(i:-1:i-M+1);
    E(i)= D(i) - (wi)*(j)' ;
    wi = wi + 2*mu*E(i)*(j);
end

% estimation of the signal
Est = zeros(n,1);
for i = M:n
   j =  A(i:-1:i-M+1);
   Est(i) = ((wi)*(j)')
end

% computing error signal
Err = Est' - D ;

subplot(4,1,1);
plot(D);
title('Desired signal');

subplot(4,1,2);
plot(A);
title('noisy signal');

subplot(4,1,3);
plot(Est);
title('Estimation signal');

subplot(4,1,4);
plot(Err);
title('Error signal');
