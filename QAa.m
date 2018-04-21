%length of input
close all;
N= 20000;
m = 1:N;
%update speed of weights
mu = 10e-8;
% not sure what this parameter is called and its function
r = 0.98;
w0=0.9*pi;

s = cos(w0*m/1000);
w_noise=0.2*pi; 
noise =100*cos(w_noise*m); 
x = s+noise;

% ADD NOISE
[y,a] = adaptive_filter(x,N,mu,r);

figure(1);
subplot(2,1,1);
stem(m,y);
xlabel('Samples(n)');
ylabel('y');
title('Output signal');

subplot(2,1,2);
stem(m,s);
xlabel('Samples(n)');
ylabel('s');
title('Desired signal');

legend('y','s');
figure(2)
plot(a);
xlabel('stages(n)');
ylabel('a');
title('Weights');
legend('a');







