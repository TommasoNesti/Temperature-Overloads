function [det,curr,lb,taylor,limits_det,limits_current,limits_lower_bound,limits_taylor]=ieee(p,eps,T,tau,gamma,L)




%load case
define_constants
%mpc=loadcase('case145');
mpc=loadcase('c118swf');
branch=mpc.branch;
bus=mpc.bus;
gen_index=mpc.gen(:,1);

%build matrices
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
b = stat ./ branch(:, BR_X);                    %% series susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
b = b ./ tap;

%% build connection matrix Cft = Cf - Ct for line and from - to buses
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
i = [(1:nl)'; (1:nl)'];                         %% double set of row indices
Cft = sparse(i, [f;t], [ones(nl, 1); -ones(nl, 1)], nl, nb);    %% connection matrix

Bf = sparse(i, [f; t], [b; -b], nl, nb);    % = spdiags(b, 0, nl, nl) * Cft;

Bbus = Cft' * Bf;

slack = find(bus(:, BUS_TYPE) == REF);%find slack bus
nb = size(bus, 1);
nbr = size(branch, 1);
noref   = (2:nb)';      %% use bus 1 for voltage angle reference
slack_bus = slack;
noslack = find((1:nb)' ~= slack_bus);

H = zeros(nbr, nb);
H(:, noslack) = full(Bf(:, noref) / Bbus(noslack, noref));

H_check=makePTDF(mpc);
norm(H-H_check);
C_tilde=H(:,2:end);


%
d=mpc.bus(:,PD);
gen_index=mpc.gen(:,1);

%solve vanilla OPF to get original line loading
mpopt = mpoption(  'out.all', 0);

results=rundcopf(mpc,mpopt);
%results.f;
%results.branch(:,RATE_A);
g_temp=results.gen(:,PG);
g=zeros(nb,1);
g(gen_index)=g_temp;
flows=results.branch(:,PF);

%K=1.0;
%limits_det=K*abs(flows);
%Line loading are the ones provided in the test case (next they will be
%scaled down)
big_number=1e10;
limits_det=results.branch(:,RATE_A);


%min(abs(flows)./limits_det)
%mean(abs(flows)./limits_det)
%median(abs(flows)./limits_det)
%max(abs(flows)./limits_det)

%%%%%%%%



%reduce line limits to spice things up

chgtab = [ 1 1 CT_TBRCH 0 RATE_A CT_REL 0.5];
mpc = apply_changes(1, mpc, chgtab);


results=rundcopf(mpc,mpopt);
%results.f;
%results.branch(:,RATE_A);
g_temp=results.gen(:,PG);
g=zeros(nb,1);
g(gen_index)=g_temp;
flows=results.branch(:,PF);

big_number=1e10;
limits_det=results.branch(:,RATE_A);

%limits_det(find(limits_det==0))=big_number;

Delta=diag(1./limits_det);
C_bar=Delta*H;




slack_index=slack; %check for 118
stoch_index=[1,6,9,18,19,41,43,62,63,72,80]';
det_index=setdiff([1:nb]',stoch_index);
nb_s=length(stoch_index);
C_bar=Delta*H;
C=C_bar(:,stoch_index);
C_D=C_bar(:,det_index);

%Initialize OU process for stochastic buses.
%Each bus is different. D,L,T,eps, sono passati da input

tau_0=tau;
%D=gamma*eye(nb_s,nb_s);
%L=l*eye(nb_s,nb_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%define new line limits according to LD expressions

%0) determinsitic region
"DETERMINISTIC"
limits_det;

mpc_det=mpc;
for line=1:nl
    f_limit=limits_det(line);
    chgtab = [ 
    1 1 CT_TBRCH line RATE_A CT_REP f_limit;
    ];
    mpc_det = apply_changes(1, mpc_det, chgtab);
end

branch_det=mpc_det.branch;
branch_det(:,RATE_A);
results_det=rundcopf(mpc_det,mpopt);
flows_det=results_det.branch(:,PF);
results_det.f



%1) Current region
"DEFINE CURRENT REGION"



%sigma_sq=diag(C*L*L*C'); %
%PROBLEMl IN THE PAPER THE SLACK BUS IS NOT COUNTED HERE
%THIS IS THE MAIN DRIVER OF HOW SMALL BETA IS. INVESTIAGE


%the old version is correct, since the only thing that changed is L
% while D is still a multiple of the identiy
%L=l*eye(nb_s,nb_s);
sigma_sq=diag(C*L*L*C'); 
num_c=eps*log(1/p)*(1-exp(-2*gamma*T));
den_c=gamma;
sqrt(num_c/den_c)
max(sqrt(sigma_sq))
beta=sqrt(num_c/den_c) * sqrt(sigma_sq);
limits_current=(1-beta).*limits_det;

%current decay rate



%new version (stochastic buses not uniform)
%D=gamma*eye(nb_s,nb_s);
%I=eye(nb_s);
%M=L*L*inv(D)*(I-expm(-2*T*D));
%sigma_sq_combined=diag(C*M*C');
%sigma_sq=diag(C*L*L*C');
%beta=sqrt(eps*log(1/p))*sqrt(sigma_sq_combined);
%limits_current=(1-beta).*limits_det;


if ~all(limits_current>0)
    "ERROR"
end

%WATCH OUT: i line limtis glieli devo passare non scalati. forse
%moltiplicare per line limits

mpc_current=mpc;
for line=1:nl
    f_limit=limits_current(line);
    chgtab = [ 
    1 1 CT_TBRCH line RATE_A CT_REP f_limit;
    ];
    mpc_current = apply_changes(1, mpc_current, chgtab);
end



branch_current=mpc_current.branch;
branch_current(:,RATE_A)';
results_current=rundcopf(mpc_current,mpopt);
flows_current=results_current.branch(:,PF);
results_current.f;


%D=gamma*eye(nb_s,nb_s);
%I=eye(nb_s);
%M=L*L*inv(D)*(I-expm(-2*T*D));
%sigma_sq_combined=diag(C*M*C');
%nu=flows_current./limits_det
%I_c=min((1-abs(nu)).^2./sigma_sq_combined)


%figure;
%plot(limits_det)
%hold on
%plot(limits_current)
%hold on
%plot(abs(flows_det))

%2) lower bound for temperature
"DEFINE LOWER BOUND REGION"

%old version
num_lb=eps*log(1/p)*(1-exp(-2*gamma*T));
den_lb=gamma;
eta=sqrt(num_lb/den_lb)* sqrt(sigma_sq); %equal to beta. CHECK
all(eta==beta)
delta=sqrt(1-eta.*eta*exp(-T/tau)*(1-exp(-T/tau))) -(1-exp(-T/tau))*eta;
limits_lower_bound =real(delta.*limits_det);% to remove



%new version
%num_lb=eps*log(1/p)*(1-exp(-2*gamma*T));
%den_lb=gamma;
%eta=sqrt(num_lb/den_lb)* sqrt(sigma_sq); 
%all(eta==beta)
%delta=sqrt(1-eta.*eta*exp(-T/tau)*(1-exp(-T/tau))) -(1-exp(-T/tau))*eta;
%limits_lower_bound=real(delta.*limits_det); % da rimettere


isreal(limits_lower_bound)
%INTRODUCE SOME CHECKS AGAINST COMPELS VALUES

if ~all(limits_lower_bound>0)
    "ERROR"
end

mpc_lower_bound=mpc;
for line=1:nl
    f_limit=limits_lower_bound(line);
    chgtab = [ 
    1 1 CT_TBRCH line RATE_A CT_REP f_limit;
    ];
    mpc_lower_bound = apply_changes(1, mpc_lower_bound, chgtab);
end

branch_lower_bound=mpc_lower_bound.branch;
branch_lower_bound(:,RATE_A);
results_lower_bound=rundcopf(mpc_lower_bound,mpopt);
flows_lower_bound=results_lower_bound.branch(:,PF);


%3) taylor expansion (REAL PARTS)
"DEFINE TAYLOR REGION"
%OLD

%limits_taylor=real((1-1/sqrt(1+2*tau_0*gamma)*eta).*limits_det);

%new
limits_taylor=real((1-eta/sqrt(1+2*tau_0*gamma)).*limits_det);
isreal(limits_taylor)

if ~all(limits_taylor>0)
    "ERROR"
end

mpc_taylor=mpc;
for line=1:nl
    f_limit=limits_taylor(line);
    chgtab = [ 
    1 1 CT_TBRCH line RATE_A CT_REP f_limit;
    ];
    mpc_taylor = apply_changes(1, mpc_taylor, chgtab);
end

branch_taylor=mpc_taylor.branch;
branch_taylor(:,RATE_A);
results_taylor=rundcopf(mpc_taylor,mpopt);
flows_taylor=results_taylor.branch(:,PF);
%print different limits
limits_det;
limits_current;
limits_lower_bound;
limits_taylor;

%print differnet objective functions
det=results_det.f
curr=results_current.f
lb=results_lower_bound.f
taylor=results_taylor.f





if ~all(limits_current>0)==1
    "CURRENT WRONG"
    min(limits_current)
end
if ~all(limits_lower_bound>0)==1
    "LB WRONG"
    min(limits_lower_bound)
end

if ~all(limits_taylor>0)==1
    "TAYLOR WRONG"
    min(limits_taylor)
end

%figure;
%hold on;
%plot(flows_det)
%plot(flows_current)
%plot(flows_lower_bound)
%plot(flows_taylor)
%legend('det','cur','lb','tayl')
%figure;
%plot(limits_det)
%hold on
%plot(limits_current)
%hold on

%plot(limits_lower_bound)
%hold on
%plot(limits_taylor)
end


