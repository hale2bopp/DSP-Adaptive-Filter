%% cascade notch
m = 1:N;
N= 20000;

%input1 = cos(pi*m); 


rand_i = rand(1,N);
rand_i(rand_i>=0.5)=1;
rand_i(rand_i<0.5)=0;

input1 = rand_i;
	
noise_1 = 100*cos(0.2*pi*m); 
noise_2 = 100*cos(0.6*pi*m); 

x = input1 + noise_1 + noise_2;

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

figure(1);

stem(m,input1);
xlabel('Samples(n)');
ylabel('Input signal (cascade)');
title('Input signal (cascade)');

hold on;

stem(m,y2);
xlabel('Samples(n)');
ylabel('Output signal (cascade)');
title('Output  signal (cascade)');

hold off;

legend('input','y(cascade)');

figure(2)
plot(a1);
xlabel('Samples(n)');
ylabel('Weights');
title('Weights');
hold on;

plot(a2);
xlabel('Samples(n)');
ylabel('Weights');
title('Weights');
hold off;

legend('a1','a2');
