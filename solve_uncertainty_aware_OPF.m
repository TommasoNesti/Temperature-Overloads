function [det,curr,lb,taylor,limits_det,limits_current,limits_lower_bound,limits_taylor]=solve_uncertainty_aware_OPF(p,eps,T,tau,gamma,L)

% This function is called by the script "generate_fig2" to reproduce 
%  Figure 2 in the paper
%
% [1] "Temperature Overloads in Power Grids Under Uncertainty:  a Large 
%      Deviations Approach", Tommaso Nesti, Jayakrishnan Nair 
%      and Bert Zwart, accepted for publication in IEEE Transactions on 
%      Control of Network Systems.
%
%Author: Tommaso Nesti
%Version 1.0 (May 2019)
%
% Notes:
% 1) This function solves four distinct DC OPFs, each incorporating a 
% different capacity region in the constraints (deterministic, 
% current-based, lower-bound based and Taylor-expansion-based.
%
% 2) The code heavily relies on MATPOWER functions
% (http://www.pserc.cornell.edu/matpower/), and uses the test case 'c118swf' 
% (https://github.com/MATPOWER/matpower/blob/master/most/lib/t/c118swf.m).
% The test case differs from the standard IEEE 118 test case in two ways:
% it incorporates line limits, and it differentiates between deterministic 
% and renewable generators.


%load case
define_constants
mpc=loadcase('c118swf');
branch=mpc.branch;
bus=mpc.bus;
gen_index=mpc.gen(:,1);
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%build PTDF matrix
H=makePTDF(mpc);

%reduce line limits by 50 %
chgtab = [ 1 1 CT_TBRCH 0 RATE_A CT_REL 0.5];
mpc = apply_changes(1, mpc, chgtab);
mpopt = mpoption(    'out.all', 0);

% solve preliminary deterministic DC OPF 
results=rundcopf(mpc,mpopt);
limits_det=results.branch(:,RATE_A);

% initialize parameters described in Section II ("System Model") of [1]
Delta=diag(1./limits_det);
stoch_index=[1,6,9,18,19,41,43,62,63,72,80]';
C_bar=Delta*H;
C=C_bar(:,stoch_index);


% Solve the 4 DC OPFs:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Determinsitic region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Current region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_sq=diag(C*L*L*C'); 
num_c=eps*log(1/p)*(1-exp(-2*gamma*T));
den_c=gamma;
sqrt(num_c/den_c)
max(sqrt(sigma_sq))
beta=sqrt(num_c/den_c) * sqrt(sigma_sq);
limits_current=(1-beta).*limits_det;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Lower bound-based region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_lb=eps*log(1/p)*(1-exp(-2*gamma*T));
den_lb=gamma;
eta=sqrt(num_lb/den_lb)* sqrt(sigma_sq); 
all(eta==beta)
delta=sqrt(1-eta.*eta*exp(-T/tau)*(1-exp(-T/tau))) -(1-exp(-T/tau))*eta;
limits_lower_bound =real(delta.*limits_det);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Taylor expansion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
limits_taylor=real((1-eta/sqrt(1+2*tau*gamma)).*limits_det);
isreal(limits_taylor)
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


% Store different objective functions for output
det=results_det.f;
curr=results_current.f;
lb=results_lower_bound.f;
taylor=results_taylor.f;
end


