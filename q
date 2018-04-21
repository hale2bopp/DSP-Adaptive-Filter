%% QC

%n1 = 0:5000;

L = 10;
N = 100000;
n = 1:N;
freq = linspace(0, 0.5, L);

%{
x = zeros(L,N/L);
i=1;
for f = freq
    x(i,:) =  cos(f*2*pi*(1:L:N)); 
    i=i+1;
end
%x=x(1:N);
x=reshape(x,1,N);
%}
x = zeros(1,N);
for f = freq
    x = x+ cos(f*2*pi*n); 
end

f = linspace(0,0.5,N);


Hd = 2*(f>=0 & f<=0.15) + 1*(f>=0.3 & f<= 0.5);

for n2 = 1:N
    Hd(n2) = abs(Hd(n2))*exp(1i*f(n2)*2*pi);
end

figure(20);
freqz(Hd,1);

%hd = ifft(Hd);
X_desired = fft(x).*Hd(1:N);
conv_h = ifft(X_desired);

M = 40;
r = 0.85;
adapt_filter = zeros(1, M+1);
adapt_filter(M/2) = 1;
y = zeros(1,N);
mu = 1e-7;

for k=M+1:N
    for indx = 1:M+1
        y(k)= y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error(k) = abs(conv_h(k) - y(k));
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) + 2*mu*error(k)*x(k+1-indx);
    end
end

figure; plot(error); title('error');
figure; freqz(adapt_filter,1);title('adaptive filter magnitude');
matched_filt = abs(adapt_filter);
figure; plot(abs(matched_filt)); xlabel('frequency'); ylabel('adaptive filter'); title('magnitude of adaptive filter');

a
