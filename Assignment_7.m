%% #1a

f=1;
A=1;

t = 0:0.001:2*pi;
x_t = A*cos(2*pi*f*t);
figure,plot(t,x_t)

%% #1b
f_s = 5;
T_s = 1/f_s;
time_step = 0.001*T_s;
t = 0:time_step:2*pi;
n = 0:10;

x_t = A*cos(2*pi*f*t);
x_n = A*cos(2*pi*f*n*T_s);

figure,plot(T_s*n,x_n)

%% #1c

h = zeros(1,2.*round(T_s/time_step));
h((1:end)*time_step<T_s+time_step) = 1;

x_s = zeros(1,length(x_t));
x_s(1:(T_s/time_step):length(x_t)) = x_n;

x0 = conv(x_s,h);
x0 = x0(1:length(x_s));
hold on

%plot zero hold
subplot(1,2,1),plot((1:length(h)-1)*time_step,h)
%plot output of zero hold
subplot(1,2,2),plot(t,x0)

%% #1d
bits = 4;
maxval = 1;
minval = -1;
t_range = minval:0.001:maxval;

step = -minval/(2^(bits-1));

x_t = maxval*cos(2*pi*f*t);
k2 = round((x_t-minval)./step);
k2(k2>2^(bits)-1) = 2^(bits)-1;
x_t_q = minval+k2*step;

figure, plot(t,x_t_q,x_t)


%% #2a

num1 = poly([0.98*exp(j*.8*pi) 0.98*exp(-j*.8*pi)]);
den1 = poly([0.8*exp(j*.4*pi) 0.8*exp(-j*.4*pi)]);;;
k = 1:4;
ck = 0.95*exp(j*(0.15*pi+0.02*pi*k));
c1 = ck(1);

amp2_1 = c1*conj(c1);
num2_1 = [1/conj(c1) 1/conj(c1) 1/c1 1/c1];
den2_1 = [c1 c1 conj(c1) conj(c1)];
c2 = ck(2);
amp2_2 = c2*conj(c2);
num2_2 = [1/conj(c2) 1/conj(c2) 1/c2 1/c2];
den2_2 = [c2 c2 conj(c2) conj(c2)];
c3 = ck(3);
amp2_3 = c3*conj(c3);
num2_3 = [1/conj(c3) 1/conj(c3) 1/c3 1/c3]; 
den2_3 = [c3 c3 conj(c3) conj(c3)];
c4 = ck(4);
amp2_4 = c4*conj(c4);
num2_4 = [1/conj(c4) 1/conj(c4) 1/c4 1/c4];
den2_4 = [c4 c4 conj(c4) conj(c4)];

amp2 = amp2_1*amp2_2*amp2_3*amp2_4;
num2 = poly([num2_1 num2_2 num2_3 num2_4]);
den2 = poly([den2_1 den2_2 den2_3 den2_4]);

num = amp2*conv(num1, num2);
den = conv(den1, den2);

%fvtool(num,den,'grpdelay')


%% #2b
n = 0:300;
M = 60;
w=@omega;
x = [];
for k=0:300
    x_1 = w(k-2*M-2)*cos(.4*pi*(k-2*M-2)-pi/2);
    x_2 = w(k-M-1)*cos(0.2*pi*(k-M-1));

    x = [x,x_1+x_2];
end

figure,plot(n,x)

y = filter(num,den,x);

figure, plot(n,y)

function w = omega(n)
    M=60;
    w = 0.54-.46*cos(2*pi*n/M);
    w = w*heaviside(n-1)*(1-heaviside(n-M));
end
