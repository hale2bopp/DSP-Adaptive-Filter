%training_sequence = normrnd(1,1000);
%h = [0.3, 1,0.7,0.3,0.2];
%gaussian_noise = awgn(x,25,'measured'); % assume snr 25dB
% updates

%length of input
close all;
N= 20000;
m = 1:N;

% weights  
a = zeros(1,N);
e = zeros(1,N);
% initialize an empty array for output (and also to set the initial value
% for output
y =zeros(1,N);


%update speed of weights
mu = 10e-8;

% not sure what this parameter is called and its function
r = 0.98;


w0=pi;
s = cos(w0*n);

w_noise=0.2*pi; 
noise =100*cos(w_noise*m); 

x = s+noise;

% ADD NOISE


for n=3:N
    e(n) = x(n) + a(n)*x(n-1) +x(n-2);
    y(n) = e(n) - r*a(n)*y(n-1) - (r^2)*y(n-2);
    a(n+1) =a(n) - mu*y(n)*x(n-1);
    if (abs(a(n+1))>2)
        a(n+1) = 0;
    end
end

%a=a(1:N);

% try visualizing the outputs
figure(1);
subplot(3,1,1);
stem(m,input);
xlabel('Samples(n)');
ylabel('Input signal');
title('Input signal');

hold on;
stem(m,y);
xlabel('Samples(n)');
ylabel('Output signal');
title('Output  signal');
hold off;

legend('input','y');
subplot(3,1,2);
stem(m,x);
xlabel('Samples(n)');
ylabel('noisy input signal');
title('noisy input signal');
legend('noisy input');
% weights
%figure(3)
%{
subplot(3,1,3);
stem(m,y);
xlabel('Samples(n)');
ylabel('Output signal');
title('Output  signal');
hold off;
%}
legend('input','noisy input ','y');

figure(2)
plot(a);
xlabel('Samples(n)');
ylabel('Weights');
title('Weights');
hold on;
% now outputs
%figure(3)
w1 = linspace(0,2*pi,N);

figure(3);
subplot(4,1,1);
plot(w1,abs(fft(input)));
xlabel('Frequency(w)'); 
ylabel('Input signal');
title('Input signal FFT');
%hold on;
X=fft(x);
subplot(4,1,2);
plot(w1,abs(X));
xlabel('Frequency(w)');
ylabel('noisy input signal');
title('noisy input signal FFT');
% weights
%figure(3)
%hold on;
Y=fft(y);
subplot(4,1,3);
plot(w1,abs(Y));
xlabel('Frequency(w)');
ylabel('Output signal');
title('Output  signal FFT');
H=Y./X;
subplot(4,1,4);
plot(w1,abs( H  ));
xlabel('Frequency(w)');
ylabel('Filter');
title('Transfer function FFT');


hold off;

%% changing frequency

k = 1e-5;
N  = 20000;
x  = zeros(1,N);
for n = 1:N
    x(n) = 100*cos(2*pi*(0.5+k*n)*n);
end


%% cascade notch
N= 20000;
m = 1:N;

input = cos(pi*m); 
noise_1 = 100*cos(0.2*pi*m); 
noise_2 = 100*cos(0.6*pi*m); 

x = input + noise_1 + noise_2;

mu = 1e-6;
r = 0.98;

a1 = zeros(1,N);
e1 = zeros(1,N);
y1 = zeros(1,N);

a2 = zeros(1,N);
e2 = zeros(1,N);
y2 = zeros(1,N);


for n=3:N
    e1(n) = x(n) + a1(n)*x(n-1) + x(n-2);
    y1(n) = e1(n) - r*a1(n)*y1(n-1) - r*r*y1(n-2);
    e2(n) = y1(n) + a2(n)*y1(n-1) + y1(n-2);
    y2(n) = e2(n) - r*a2(n)*y2(n-1) - r*r*y2(n-2);
    
    a1(n+1) = a1(n) - mu*y2(n)*x(n-1);
    
    if (abs(a1(n+1))>2)
        a1(n+1) = 0;
    end
 
    a2(n+1) = a2(n) - mu*y2(n)*y1(n-1);
    
    if (abs(a2(n+1))>2)
        a2(n+1) = 0;
    end
    
end

figure; plot(a1);
figure; plot(a2);


%% Q B

N = 5000;

train = sign(rand(1, N)-0.5);
train = 20*train;

h = [0.3,1,0.7,0.3,0.2];
noise = randn( 1, length(conv(h,train)) );
x = conv(h, train) + noise;


d = 2;
M = 20;
adapt_filter = zeros(1,M+1);
adapt_filter(M/2) = 1;
mu = 1e-6;
error = zeros(1,N);
y = zeros(1,N);

for k=M+1:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = train(k-d) - y(k);
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) + 2*mu*error(k)*x(k+1-indx);
    end
end

figure; plot(error);
figure; freqz(h,1);
figure; freqz(adapt_filter,1);




%% Q B Part 2

N = 100000;

train = sign(rand(1, N)-0.5);
train = 100*train;

h = [0.3,1,0.7,0.3,0.2];

%noise = randn( 1, length(conv(h,train)))/10;
noise = randn( 1, length(conv(h,train)));
x = conv(h, train) + noise;

M = 20;
adapt_filter = zeros(1,M+1);
adapt_filter(M/2) = 1;
y = zeros(1,N);
%y(M/2) = 0.1;
error = zeros(1,N);
mu = 1e-11;

for k=M+1:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = y(k)*y(k) - 100*100 ;
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) - mu*error(k)*y(k)*x(k+1-indx);
    end
end

figure; plot(error);
figure; freqz(h,1);
figure; freqz(adapt_filter,1);

%% QC

%n1 = 0:5000;

L = 10;
N = 50000;
n = 1:N;
freq = linspace(0, 0.5, L);
x = zeros(1,N);
for f = freq
    x = x + cos(f*pi*n);
end

f = 0:0.00001:0.5;

Hd = 2*(f>=0 & f<=0.15) + 1*(f>=0.3 & f<= 0.5);

for n2 = 1:N
    Hd(n2) = abs(Hd(n2))*exp(1i*n2*pi/20);
end

hd = ifft(Hd);
conv_h = conv(hd, x);

M = 40;
N = 50000;
adapt_filter = zeros(1, M+1);
adapt_filter(M/2) = 1;
y = zeros(1,N);
mu = 1e-7;

for k=M+1:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = abs(conv_h(k) - y(k));
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) + 2*mu*error(k)*x(k+1-indx);
    end
end

figure; plot(error);
figure; freqz(adapt_filter,1);

%% QC part 2


L = 10;

freq = linspace(0, 0.5, L);
x = zeros(1,N);
N = 50000;
n = 0:N;
for f = freq
    x = x + 100*cos((pi+f)*n);
end

f = 0:0.0001: 0.5;

Hd = 2*pi*f*(f>=0.0 & f<=0.30);

for n2 = 1:N
    Hd(n2) = abs(Hd(n2))*exp(1i*(n2*pi/20 + pi/2));
end

hd = ifft(Hd);
conv_h = conv(hd, x);

M = 40;
N = 50000;
adapt_filter = zeros(1, M+1);
adapt_filter(M/2) = 1;
y = zeros(1,N);
mu = 1e-7;

for k=M+1:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = abs(conv_h(k) - y(k));
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) + 2*mu*error(k)*x(k+1-indx);
    end
end









