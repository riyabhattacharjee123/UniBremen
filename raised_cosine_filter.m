close all;
clear all;
clc;

L=41; %filter Length
Data_rate = 72.25; % ; % bps Y1 ; %data Rate = 2.5Gbps
Fs=8*Data_rate; %oversampling by 8
T=1/Data_rate; %pulse duration
Ts=1/Fs; %sampling duration
alpha =1; % Design Factor for Raised Cosing Filter
upsampling = Fs/Data_rate;


%data1 = round(rand(64,2^3)); % Generate a random bit stream
bits = 2^5;

%----------------------------------------------------------
%Raised Cosing Filter Design
%----------------------------------------------------------
if mod(L,2)==0 
 M=L/2 ; % for even value of L
else
 M=(L-1)/2; % for odd value of L 
end 
g=zeros(1,L); %Place holder for RC filter's transfer function 
for n=-M:M 
 num=sin(pi*n*Ts/T)*cos(alpha*pi*n*Ts/T);
 den=(pi*n*Ts/T)*(1-(2*alpha*n*Ts/T)^2);
 g(n+M+1)=num/den; 
 if (1-(2*alpha*n*Ts/T)^2)==0 
 g(n+M+1)=pi/4*sin(pi*n*Ts/T)/(pi*n*Ts/T);
 end 
 if (pi*n*Ts/T)==0 
 g(n+M+1)=cos(alpha*pi*n*Ts/T)/(1-(2*alpha*n*Ts/T)^2); 
 end 
end

%y11=zeros[][];
for d = 1:64
    data1 = randi([0,1],d,bits); % Generate a random bit stream
    output1(d,:)=upsample(data1(d,:),upsampling);  
end

for d=1:64
    y11(d,:)=conv(g,output1(d,:));
end



%Plot data and RC filtered Output
figure;
subplot(3,1,1);
stem(data1(1,:));
title('Input data to the Raised Cosine Filter by first GS antenna');
xlabel('Samples');
ylabel('Amplitude');

subplot(3,1,2);
stem(output1(1,:));
title('Oversampled Input data to the Raised Cosine Filter by first GS antenna');
xlabel('Samples');
ylabel('Amplitude');

subplot(3,1,3);
plot(y11(1,:));
title('Response of the Raised Cosine Filter for the given Input by first GS antenna');
xlabel('Samples');
ylabel('Amplitude');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = y11 ; % desired signal
n_size = size(D,2);
for ns = 1:3
    A = D(ns,1:n_size) + 0.9*randn(ns,n_size);
end

mmse_length = 64 ; %length_MMSE_adapter
Weight = zeros(3,mmse_length);
wi = zeros(3,mmse_length);
E = [];
step_size = 0.0005;


for idx = mmse_length:n_size
    ix = idx:-1:idx-mmse_length+1 ;
    mult(idx,:) = (wi(:,idx))'.*A
    E(idx,idx) = D(idx,idx) - mult  
    
end

