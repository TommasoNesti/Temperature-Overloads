
cd /export/scratch1/nesti/MATLAB/matpower7.0b1
startup

%RANGE OF OVERLOAD PROBABILITY
p_min=1e-1;
p_max=1e-7;
p_num=7;
exp_min=log10(p_min);
exp_max=log10(p_max);
range2=10.^(linspace(exp_min,exp_max,p_num));


%%%%%%%
base=[65,626,742,212,63,512,66,153,87,45,40]';
nb_s=length(base);
gamma=1; %0.1
tau=0.5;
eps=1;


m=3;
n=size(range2,2);
perc_vec=zeros(nb_s,3);
perc_vec=[0.01,0.018,0.04]; %CAMBIARE
var_vec=zeros(nb_s,3);

vec_det=zeros(m,n);
vec_curr=zeros(m,n);
vec_lb=zeros(m,n);
vec_taylor=zeros(m,n);
i=1;

range_T=[1/12, 1/4, 1]
for T = range_T
    j=1;

    perc=perc_vec(i);
    var=(perc*base).^2;
    var_vec(:,i)=var;
    l=sqrt(2*var*gamma/(eps*(1-exp(-2*gamma*T))));
    L=diag(l);

   for p=range2
        [det,curr,lb,taylor,limits_det,limits_current,limits_lower_bound,limits_taylor]=ieee_not_uniform(p,eps,T,tau,gamma,L);
        if length(det)==0
            det=nan;
        end
    
        if length(curr)==0
          curr=nan;
        end
    
        if length(lb)==0
            lb=nan;
        end
    
        if length(taylor)==0
            taylor=nan;
        end
        vec_det(i,j)=det;
        vec_curr(i,j)=curr;
        vec_lb(i,j)=lb;
        vec_taylor(i,j)=taylor;
        j=j+1;
       end
    i=i+1;
    
end
    
        
    %vec_curr
    %vec_lb
    
    
 figure;
 plot(limits_current)
 hold on;
 plot(limits_lower_bound)
 
 
range1=range_T;

% [X,Y] = meshgrid(range1,range2);
% X=60*X
% %X=100*X
% figure;
% %ax1 = subplot(1,1,1); 
% Z_curr=(vec_curr-vec_det)./vec_det;
% Z_curr=1*Z_curr';
% surf(X,Y,Z_curr,'FaceColor',[0 0 1],'FaceAlpha',0.5)
% set(gca,'YScale','log')

%s = surf(X,Y,Z,'FaceAlpha',0.5)

%plot3(X,Y,Z,'r')
%stem3(X,Y,Z,'o')
%colormap(ax1,winter)
% hold on;
% 
% %ax2 = subplot(1,1,1); 
% Z_lb=(vec_lb-vec_det)./vec_det;
% Z_lb=1*Z_lb';
% surf(X,Y,Z_lb,'FaceColor',[0 1 1],'FaceAlpha',0.5)
% %s = surf(X,Y,Z,'FaceAlpha',0.5)%plot3(X,Y,Z,'r')
% %stem3(X,Y,Z,'b')
% %colormap(ax2,summer)
% hold on;

%ax3 = subplot(1,1,1); 
% Z_taylor=(vec_taylor-vec_det)./vec_det;
% Z_taylor=1*Z_taylor';
% surf(X,Y,Z_taylor,'FaceColor',[1 0 0],'FaceAlpha',0.5)
% %s = surf(X,Y,Z,'FaceAlpha',0.5)%plot3(X,Y,Z,'r')
% %stem3(X,Y,Z,'r')
% %colormap(ax3,spring)
% 
% 
% set(gca,'YScale','log')
% %xlabel('OU long term Standard Deviation (% of nominal value)')
% xlabel('Time interval (minutes)')
% ylabel('Overload probability')
% 
% %zlabel({'Relative increase in system cost (%)';'(wrt. deterministic system'})
% zlabel('CoU')
% %colormap(spring)
% %legend(['Current-based' newline ' capacity region'],'Lower bound-based capacity region','Taylor-based capacity  region')
% legend('Current','Lower bound','Taylor')
% 
% %saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/CoU_118_new.png')
% 
% 
% 
% %benchmark is current
% 
% Z_lb=(vec_lb-vec_curr)./vec_curr;
% Z_lb=abs(Z_lb');
% Z_taylor=(vec_taylor-vec_curr)./vec_curr;
% Z_taylor=abs(Z_taylor');
% 
% %savings in costs
% Z_lb=100*Z_lb;
% Z_taylor=100*Z_taylor;
% 
% figure;
% semilogx(range2',Z_lb(:,1),':k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range2',Z_taylor(:,1),'--k.','LineWidth',2,'MarkerSize',20);
% xlabel('Overload probability p')
% ylabel('Reduction in system costs (%)')
% legend('Lower bound','Taylor')
% %saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/isaving_118_T=5min.png')
% 
% 
% figure;
% semilogx(range2',Z_lb(:,2),':k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range2',Z_taylor(:,2),'--k.','LineWidth',2,'MarkerSize',20);
% xlabel('Overload probability p')
% ylabel('Reduction in system costs (%)')
% legend('Lower bound','Taylor')
% %saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/saving_118_T=5min.png')
% 
% figure;
% semilogx(range2',Z_lb(:,3),':k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range2',Z_taylor(:,3),'--k.','LineWidth',2,'MarkerSize',20);
% xlabel('Overload probability p')
% ylabel('Reduction in system costs (%)')
% legend('Lower bound','Taylor')
% %saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/saving_118_T=60min.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%benchmark is deterministic
Z_curr=(vec_curr-vec_det)./vec_det;
Z_curr=Z_curr';
Z_lb=(vec_lb-vec_det)./vec_det;
Z_lb=Z_lb';
Z_taylor=(vec_taylor-vec_det)./vec_det;
Z_taylor=Z_taylor';

%Z_curr=100*Z_curr;
%Z_lb=100*Z_lb;
%Z_taylor=100*Z_taylor;

% figure;
% semilogx(range2',Z_curr(:,1),'-k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range2',Z_lb(:,1),':k.','LineWidth',2,'MarkerSize',20);
% hold on;
% semilogx(range2',Z_taylor(:,1),'--k.','LineWidth',2,'MarkerSize',20);
% xlabel('Overload probability p')
% %ylabel('Increase in system costs (%)')
% ylabel('CoU')
% legend('Current','Lower bound','Taylor')
%saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/increase_118_T=5min.png')


figure;
semilogx(range2',Z_curr(:,2),'-k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range2',Z_lb(:,2),':k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range2',Z_taylor(:,2),'--k.','LineWidth',2,'MarkerSize',20);
xlabel('Overload probability p')
%ylabel('Increase in system costs (%)')
ylabel('CoU')

legend('Current','Lower bound','Taylor')
%saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/increase_118_T=15min.png')



figure;
semilogx(range2',Z_curr(:,3),'-k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range2',Z_lb(:,3),':k.','LineWidth',2,'MarkerSize',20);
hold on;
semilogx(range2',Z_taylor(:,3),'--k.','LineWidth',2,'MarkerSize',20);
xlabel('Overload probability p')
%ylabel('Increase in system costs (%)')
ylabel('CoU')

legend('Current','Lower bound','Taylor')
%saveas(gcf,'/export/scratch1/nesti/Dropbox/Work/CWI/Resubmission_NestiZwartNair/ResubmissionIEEEcones_2018/increase_118_T=60min.png')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure;
semilogx(range2',Z_curr(:,1),'--^k','LineWidth',2);
hold on;
semilogx(range2',Z_lb(:,1),'--ok','LineWidth',2);
hold on;
semilogx(range2',Z_taylor(:,1),'--dk','LineWidth',2);
xlabel('Overload probability')
ylabel('CoU')
%legend('Current','Lower bound','Taylor')


%figure;
semilogx(range2',Z_curr(:,2),':^k','LineWidth',2);
hold on;
semilogx(range2',Z_lb(:,2),':ok','LineWidth',2);
hold on;
semilogx(range2',Z_taylor(:,2),':Xk','LineWidth',2);
xlabel('Overload probability')
ylabel('CoU')
%legend('Current','Lower bound','Taylor')




semilogx(range2',Z_curr(:,3),'-^k','LineWidth',2);
hold on;
semilogx(range2',Z_lb(:,3),'-ok','LineWidth',2);
hold on;
semilogx(range2',Z_taylor(:,3),'-Xk','LineWidth',2);
xlabel('Overload probability')
ylabel('CoU')

legend('Current','Lower bound','Taylor')


Z_curr
Z_lb
Z_taylor