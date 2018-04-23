%% QC

%n1 = 0:5000;

L = 10;
N = 100000;
n = 1:N;
freq = linspace(0, 0.5, L);

t = n./1000;
x = 20*zeros(1,N);
for f = freq
    x = x+ cos(f*2*pi*t); 
end

f = linspace(0,0.5,N);


Hd = 2*(f>=0 & f<=0.15) + 1*(f>=0.3 & f<= 0.5);

for n2 = 1:N
    Hd(n2) = abs(Hd(n2))*exp(-1*1i*f(n2)*2*pi);
end


hd = ifft(Hd);
X_desired = fft(x).*Hd(1:N);
conv_h = ifft(X_desired);

%conv_h = conv(hd,x);
M = 100;
r = 0.98;
adapt_filter = zeros(1, M+1);
adapt_filter(M/2) = 1;
y = zeros(1,length(x)+N+1);
mu = 1e-7;

for k=M+1:N
    for indx = 1:M+1
        y(k)= y(k) + x(k+1-indx)*adapt_filter(indx);
    end
    error1(k) = conv_h(k) - y(k);
    for indx = 1:M+1
        adapt_filter(indx) = adapt_filter(indx) + 2*mu*error1(k)*x(k+1-indx);
    end
end

freq_vect = linspace(0,0.5,length(adapt_filter));
figure; plot(abs(error1)); title('error');
figure; freqz(adapt_filter,1);title('adaptive filter magnitude');
matched_filt = abs(adapt_filter);
figure; plot(freq_vect,matched_filt); xlabel('frequency'); ylabel('adaptive filter'); title('magnitude of adaptive filter and desired filter');hold on;
plot(f,abs(Hd)); xlabel('frequency'); ylabel('Desired filter');hold off;
%figure; subplotstem(y);xlabel();


