%% changing frequency

k = 1e-5;
N  = 20000;
x  = zeros(1,N);

mu = 10e-8;
r = 0.98;

for n = 1:N
    x(n) = 100*cos(2*pi*(0.5+k*n)*n);
end

[y,a] = adaptive_filter(x,N,mu,r);

figure(1);
subplot(2,1,1);
stem(m,y);
xlabel('Samples(n)');
ylabel('y');
title('Output signal');

subplot(2,1,2);
stem(m,x);
xlabel('Samples(n)');
ylabel('s');
title('Desired signal');

legend('y','s');

plot(a);
xlabel('Samples(n)');
ylabel('noisy input signal');
title('noisy input signal (cascade)');
legend('noisy input');


