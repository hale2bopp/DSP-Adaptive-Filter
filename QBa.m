%% Q B
N = 5000;

train = sign(rand(1, N)-0.5);
train = 20*train;

h = [0.4,0.5,1,0.2,0.3];
noise = randn( 1, length(conv(h,train)) );
x = conv(h, train) + noise;

d = 2;

M = 20;
adapt_filter = zeros(1,M+1);
adapt_filter(ceil(M/2)) = 1;
mu = 1e-6;
error = zeros(1,N);
y = zeros(1,N);

for k=M+1:N
    for indx = 1:M+1
        y(k) = y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = (train(k-d)) - y(k);
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) + 2*mu*error(k)*x(k+1-indx);
    end
end

figure; plot(error); title('Error h = [0.4,0.5,1,0.2,0.3]');
figure; freqz(h,1);title('Channel Impulse Response h = [0.4,0.5,1,0.2,0.3]');
figure;freqz(adapt_filter,1);title('Adaptive Filter Frequency Response h = [0.4,0.5,1,0.2,0.3]');


% conv(h,adapt_filter) should have a notch response

