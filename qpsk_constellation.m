 clear;
close all;
nr_data_bits=64000;
nr_symbols=nr_data_bits/2;
b_data=(randn(1,nr_data_bits)>.5);
b=[b_data];

d=zeros(1,length(b)/2);

for n=1:length(b)/2;
    p=b(2*n);
    imp=b(2*n-1);
    if (imp==0)&(p==0)
        d(n)=exp(j*pi/4);
    end
    if(imp==1)&(p==0)
        d(n)=exp(j*3*pi/4);
    end
    if (imp==1)&(p==1)
        d(n)=exp(j*5*pi/4);
    end
    if (imp==0)&(p==1)
        d(n)=exp(j*7*pi/4);
    end
end
qpsk=d;
figure(1);
plot(d,'o');
axis([-2 2 -2 2]);
grid on;
xlabel('real');
ylabel('imag');
title('QPSK constellation')