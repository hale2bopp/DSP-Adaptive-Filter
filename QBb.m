%% Q B Part 2
A = 20;
N = 100000;
train = sign(rand(1, N)-0.5);
train = A*train;

h = [0.3,1,0.7,0.3,0.2];

%noise = randn( 1, length(conv(h,train)))/10;
noise = randn( 1, length(conv(h,train)));
x = conv(h, train) + noise;

M = 20;
adapt_filter = zeros(1,M+1);
adapt_filter(M/2) = 1;
y = zeros(1,N);
%y(M/2) = 0.1;
error1 = zeros(1,N);
mu = 1e-10;

for k=M+1:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error1(k) = y(k)*y(k) - A*A;
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) - mu*error1(k)*y(k)*x(k+1-indx);
    end
end

figure; plot(error1);xlabel('iterations');ylabel('error');title('blind equalization error');
figure; freqz(h,1);xlabel('iterations');ylabel('error');title('blind equalization error');
figure; freqz(adapt_filter,1);title('blind equalization adaptive filter');

%convolve the adaptive filter with the channel impulse
overall_filter = conv(h,adapt_filter);
figure; freqz(overall_filter,1);title('Combined frequency response');
figure; stem(y);hold on; stem(train);
xlabel('samples');ylabel('output');title('Output sequence');legend('output','train'); 











