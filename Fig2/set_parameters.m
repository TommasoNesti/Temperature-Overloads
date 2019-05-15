%run this
%cd /export/scratch1/nesti/MATLAB/matpower7.0b1
%startup

eps=0.5;
p=1e-7;
tau=1;
T=1;
l=45;
gamma=1;
K=1.5;

var=200
l=sqrt(2*var*gamma/eps)




base=[65,626,742,212,63,512,66,153,87,45,40]';
nb_s=length(base);
perc=0.05
var=(perc*base).^2
gamma=1;
%l=1.3
T=1
eps=1;
%I allow for variable volatility l, but in the range [1,2]
l=sqrt(2*var*gamma/(eps*(1-exp(-2*gamma*T))));
tau=0.5;
p=0.0001
L=diag(l).*eye(nb_s);
%sqrt(eps*l.^2*(1-exp(-2*T*gamma))/(2*gamma))./base
1-exp(-T/tau)


%%%%%%%%
%%%%%%%
base=[65,626,742,212,63,512,66,153,87,45,40]';
nb_s=length(base);
gamma=0.13;
eps=1;
T=0.1%1/12;
perc=0.01;
tau=0.1

var=(perc*base).^2;
%var=(perc*mean(base)).^2;

l=sqrt(2*var*gamma/(eps*(1-exp(-2*gamma*T))));
sqrt(eps*l.^2*(1-exp(-2*T*gamma))/(2*gamma))
L=diag(l).*eye(nb_s); %riguardare

%sqrt(eps*l^2*(1-exp(-2*T*gamma))/(2*gamma))

%provare a fissare l e calibrare T

p=0.0001

1-exp(-T/tau)

figure;
plot(limits_lower_bound-limits_current)

figure;
plot(delta)
hold on
plot(1-beta)
