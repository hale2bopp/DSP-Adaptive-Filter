%training_sequence = normrnd(1,1000);
%h = [0.3, 1,0.7,0.3,0.2];
%gaussian_noise = awgn(x,25,'measured'); % assume snr 25dB
% updates

%length of input
close all;
N= 2000;

% weights  
a = zeros(1,N);

%update speed of weights
mu = 0.001;

% not sure what this parameter is called and its function
r = 0.85;

%initializing the error function
e = zeros(1,N);


Tfinal= 1;%s
%m = linspace(0,Tfinal,N);
m = 1:N;
%f= 105;%Hz
w=0.1*pi;
%input 
input = 2*sin(w*m);
f_noise = 875;
w_noise=0.9*pi;
%noise 
noise =3*sin(w_noise*m); % normrnd(0,1)/2; %divided by 2 

% initialize an empty array for output (and also to set the initial value
% for output
y =zeros(1,N);


% ADD NOISE

x=input+noise;

a(1) =0;
a(2) = a(1);
%a(3) = a(2) - mu*y(2)*x(1);
e(1) = x(1);
y(1) = e(1);
e(2) = x(2)+a(2)*x(1);
y(2) = e(2)-r*a(2)*y(1);
a(3) = a(2) - mu*y(2)*x(1);

%num_iter = 1000;

for n=3:N
    e(n) = x(n) + a(n)*x(n-1) +x(n-2);
    y(n) = e(n) - r*a(n)*y(n-1) - (r^2)*y(n-2);
    a(n+1) =a(n) - mu*y(n)*x(n-1);
    if ((a(n)>=2)||(a(n)<-2))
        a(n) = 0;
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

