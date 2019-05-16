% This script reproduces Figure 2 in the paper
%
% [1] "Temperature Overloads in Power Grids Under Uncertainty:  a Large 
%      Deviations Approach", Tommaso Nesti, Jayakrishnan Nair 
%      and Bert Zwart, accepted for publication in IEEE Transactions on 
%      Control of Network Systems.

%Author: Tommaso Nesti
%Version 1.0 (May 2029)

% Notes:
% 1) The values for standard deviations at the end time T of the OU process
%    are taken from:
%    [2] "B. Hodge and M. Milligan, “Wind power forecasting error 
%          distributions over multiple timescales,” in 2011 IEEE Power 
%          and Energy Society General Meeting, July 2011."
%
% 2) This script requires MATPOWER (http://www.pserc.cornell.edu/matpower/) 
%
% 3) The values for the modified 118 network, including the installed capacities, are based on test case 
%    "c118swf.m" (https://github.com/MATPOWER/matpower/blob/master/most/lib/t/c118swf.m)
%
%

cd /export/scratch1/nesti/MATLAB/matpower7.0b1 %modify with your MATPOWER folder
startup

% Range of overload probabilities
p_min=1e-1;
p_max=1e-7;
p_num=7;
exp_min=log10(p_min);
exp_max=log10(p_max);
range_p=10.^(linspace(exp_min,exp_max,p_num));

% Range of time intervals 
range_T=[1/12, 1/4, 1]; 
m=3;
n=size(range_p,2);

% standard deviations for three different time intervals: 5,15 and 60
% minutes, based on [2]
perc_vec=[0.01,0.018,0.04]; 

% installed capacities
base=[65,626,742,212,63,512,66,153,87,45,40]'; 
nb_s=length(base);

% parameters
gamma=1; 
tau=0.5;
eps=1;

%initialization
var_vec=zeros(nb_s,3);
vec_det=zeros(m,n);
vec_curr=zeros(m,n);
vec_lb=zeros(m,n);
vec_taylor=zeros(m,n);
i=1;

for T = range_T %loop over time intervals
    j=1;
    perc=perc_vec(i);
    var=(perc*base).^2;
    var_vec(:,i)=var;
    l=sqrt(2*var*gamma/(eps*(1-exp(-2*gamma*T)))); % solve Eq. (V.2) in the paper for l
    L=diag(l);

   for p=range_p %loop over overloads probabilities
       
       %solve four distinct DC OPFs, each incorporating a different
       %capacity region in the constraints (deterministic, current-based, 
       %lower-bound based and Taylor-expansion-based.
       [det,curr,lb,taylor,limits_det,limits_current,limits_lower_bound,limits_taylor]=solve_uncertainty_aware_OPF(p,eps,T,tau,gamma,L);
       vec_det(i,j)=det;
       vec_curr(i,j)=curr;
       vec_lb(i,j)=lb;
       vec_taylor(i,j)=taylor;
       j=j+1;
       end
    i=i+1;
end
    
range1=range_T;

% compute CoU, as defined in Eq. (V.3)
Z_curr=(vec_curr-vec_det)./vec_det;
Z_curr=Z_curr';
Z_lb=(vec_lb-vec_det)./vec_det;
Z_lb=Z_lb';
Z_taylor=(vec_taylor-vec_det)./vec_det;
Z_taylor=Z_taylor';


p_min_plot=0;
p_max_plot=1e-8;
p_num_plot=9;
exp_min_plot=log10(p_min);
exp_max_plot=log10(p_max);
range_p_plot=10.^(linspace(exp_min_plot,exp_max_plot,p_num_plot));


% (T=5 minutes)
% figure;
% semilogx(range_p',Z_curr(:,1),'-k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range_p',Z_lb(:,1),':k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range_p',Z_taylor(:,1),'--k.','LineWidth',2,'MarkerSize',20);
% xlabel('Overload probability p')
% ylabel('CoU')
% legend('Current','Lower bound','Taylor')

%Fig 2 (a) (T=15 minutes)
figure;
semilogx(range_p',Z_curr(:,2),'-k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range_p',Z_lb(:,2),':k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range_p',Z_taylor(:,2),'--k.','LineWidth',2,'MarkerSize',20);
%set(gca,'xtick',range_p_plot)
xlabel('Overload probability p')
ylabel('CoU')
legend('Current','Lower bound','Taylor')


% Fig 2 (b) (T=60 minutes)
figure;
semilogx(range_p',Z_curr(:,3),'-k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range_p',Z_lb(:,3),':k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range_p',Z_taylor(:,3),'--k.','LineWidth',2,'MarkerSize',20);
xlabel('Overload probability p')
ylabel('CoU')
legend('Current','Lower bound','Taylor')




