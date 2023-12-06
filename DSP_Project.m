clear all
% part 1
m=2;
f_s = 550000;
u = 3.6*10^-6;
phi_s = 0;

A_s = 1;
t = linspace(0,1*10^-3,100000);
s_t = A_s.*(t.^m.*exp(-t./u).*cos((2*pi*f_s).*t+phi_s)).*(t>0);

figure,plot(t,s_t)
xlabel('time')
xlim([0 1*10^-4])
ylabel('s_t')

h = 0.174;
SoS = 1390;
Tf= 2*h/SoS;
A_r = 0.8;  

u_r = ones(size(t));
u_r(t<Tf)=0;

r_t =s_t + (A_r/A_s*((t-Tf).^m.*exp(-(t-Tf)/u).*cos((2*pi*f_s)*(t-Tf)+phi_s)).*u_r);



figure,plot(t,r_t)
xlim([0 .5*10^-3])
xlabel('time')
ylabel('r_t')

f = 4*10^6;
T=1/f;
n = 1:1:2000;

u_r_n = ones(size(n));
u_r_n((n*T)<Tf) = 0;

s_n = A_s*(n*T).^m.*exp(-(n*T/u)).*cos((2*pi*f_s)*(n*T)+phi_s).*(n*T>0);

rcv_n = s_n + A_r/A_s*(n*T-Tf).^m.*exp(-(n*T-Tf)/u).*cos(2*pi*f_s*(n*T-Tf)+phi_s).*u_r_n;

rcv_e = A_s*(n).^m.*exp(-n/u)+A_r*(n-Tf/T).^m.*exp(-(n-Tf/T)/u).*u_r_n;

figure,plot(n,rcv_n)
xlabel('samples')
ylabel('rcv_e')

ToF = 1009*T*SoS/2;

%% part 2
noise = 0.4*cos(2*pi*f_s*n*T);
rcv_noise = rcv_n + noise.*(n*T>0);

rcv_hilbert = hilbert(rcv_noise);

figure,plot(n,abs(rcv_hilbert))
xlabel('samples')
ylabel('rcv_hilbert')

difference = 1031 - 29;

ToF_2 = difference*T*SoS/2



%% part 3
k = 1001;
f = 4*10^6;
lambda = -0.5:0.1:5;
error = zeros(length(lambda));
dis_error = zeros(length(lambda));

for i = 1:length(lambda)
    l = lambda(i);
    Tf = (k+l)*T;
    Tf_hat = k*T;
    error(i) = Tf_hat-Tf;
    dis = Tf*SoS/2;
    dis_error(i) = error(i)*SoS/2;
end
figure,plot(lambda,dis_error)
xlabel('lambda')
ylabel('distance error')





