clear all
close all
clc


data=csvread("states_data.csv",1,0);
err_record=[];
parr_record=[];
Deceased_record=[];
Infected_record=[];
observer=0;
%SIQRD states
S=data(:,1);
I=data(:,2);
Q=data(:,3);
R=data(:,4);
D=data(:,5);
pol=data(3:342,6:6+10-1);


par_len=9;

N=mean(S+I+Q+R+D);

l_h=1;
h_h=342;

Sp=S(l_h:h_h);
Ip=I(l_h:h_h);
Qp=Q(l_h:h_h);
Rp=R(l_h:h_h);
Dp=D(l_h:h_h);


observ_len=3;

for observe=1:h_h-observ_len+1


S=Sp(observe:observe+observ_len-1);
I=Ip(observe:observe+observ_len-1);
Q=Qp(observe:observe+observ_len-1);
R=Rp(observe:observe+observ_len-1);
D=Dp(observe:observe+observ_len-1);
days=[0 1 2];




Fs=20; 
t=linspace(0,length(days)-1,length(days)*Fs);

S_e=spline(days,S,t);
I_e=spline(days,I,t);
Q_e=spline(days,Q,t);
R_e=spline(days,R,t);
D_e=spline(days,D,t);

if (min(S_e)<0)
    S_e=S_e+abs(min(S_e));
end
if (min(I_e)<0)
    I_e=I_e+abs(min(I_e));
end
if (min(Q_e)<0)
    Q_e=Q_e+abs(min(Q_e));
end
if (min(R_e)<0)
    R_e=R_e+abs(min(R_e));
end
if (min(D_e)<0)
    D_e=D_e+abs(min(D_e));
end


S_e_d=diff(S_e)./diff(t);
I_e_d=diff(I_e)./diff(t);
Q_e_d=diff(Q_e)./diff(t);
R_e_d=diff(R_e)./diff(t);
D_e_d=diff(D_e)./diff(t);





S_e=S_e(1:end-1);
I_e=I_e(1:end-1);
Q_e=Q_e(1:end-1);
R_e=R_e(1:end-1);
D_e=D_e(1:end-1);
t=t(1:end-1);


lamb1=0.0;
lamb2=0.0000;



tries=100;
P_record=zeros(tries,par_len);
F_record=zeros(tries,1);

minner=-10;
maxxer=10;

for iter=1:tries

    p0=[-minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        0+1*rand
        ];

keta=2;
C_aditioner=1;


% zeta=p1
% beta1=p2
% beta2=p3
% alpha1=p4
% alpha2=p5
% delta=p6
% K=p7
% to=p8
% gamma=p9

N=mean(S_e+I_e+Q_e+R_e+D_e);

obj=@(p)(sum((abs(S_e_d-C_aditioner*(-p(1)*(N-S_e)+p(9)*R_e))).^keta)+...
        sum((abs(I_e_d-C_aditioner*(p(1)*(N-S_e)-p(2)*(p(8)-t).*I_e-p(3)*I_e-p(7)-p(6)*I_e))).^keta)+...
        sum((abs(Q_e_d-C_aditioner*(p(3)*I_e+p(7)-p(5)*Q_e-p(4)*Q_e))).^keta)+...
        sum((abs(R_e_d-C_aditioner*(p(2)*(p(8)-t).*I_e+p(4)*Q_e-p(9)*R_e))).^keta)+...
        sum((abs(D_e_d-C_aditioner*(p(6)*I_e.*D_e+p(5)*log(1+Q_e)))).^keta))/N;

A = [];
b = [];
Aeq = [];
beq = [];

l_lub=-50;
h_hub=50;
lb = [l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,0];
ub = [h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,1];
options = optimoptions(@fmincon, 'Algorithm' , 'interior-point','Display','off' );
options.Algorithm = 'sqp' ;
% opts.Algorithm = 'dual-simplex' ;
opts.Algorithm = 'interior-point-legacy' ;

[p,fval]=fmincon(obj,p0,A,b,Aeq,beq,lb,ub,[],options);

P_record(iter,:)=p;
F_record(iter)=fval;
end

[M,ind] = min(F_record);
p_estim=P_record(ind,1:par_len);

parr_record=[parr_record;p_estim];


omegaker=0;
if min(D_e_d)<0 && abs(min(D_e_d))>abs(max(D_e_d))
    omegaker=1;
end

omegaker2=0;
if min(I_e_d)<0 && abs(min(I_e_d))>abs(max(I_e_d))
    omegaker2=1;
end

Deceased_record=[Deceased_record;omegaker];
Infected_record=[Infected_record;omegaker2];


disp("Error:")
disp(M)

end

data_f=[parr_record Deceased_record Infected_record];
namer="new_siqrd_causalities.csv";

writematrix(data_f,namer);


