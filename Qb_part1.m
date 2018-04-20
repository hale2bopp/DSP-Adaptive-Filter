%% Q B

x = sign(rand(1, 1000)-0.5);
h = [0.3,1,0.7,0.3,0.2];
x = awgn(conv(x,h), 25);

d = 2;
M = 20;
N = 20000;
adapt_filter = zeros(1,M+1);
% what is mu?
mu = 1e-7;
error = zeros(1,N);
y = zeros(1,N);
for k=21:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = w(k-d) - y(k);
    for indx = 1:M+1
        adapt_filter(indx+1) = adapt_filter(indx) + 2*mu*error(k)*x(k+1-indx);
    end
end

