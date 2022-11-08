clear all
close all
clc


data=csvread("states_data.csv",1,0);
err_record=[];
parr_record=[];
observer=0;
%SIQRD states
S=data(:,1);
I=data(:,2);
Q=data(:,3);
R=data(:,4);
D=data(:,5);
pol=data(3:342,7:6+10);


par_len=11;

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



tries=50;
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
        -minner+maxxer*rand
        -minner+maxxer*rand
        -minner+maxxer*rand
        ];
%     p0 = (200)*rand(par_len,1);
keta=2;
E=1;
C_aditioner=0;

obj=@(p)(sum((abs(S_e_d-C_aditioner*(-p(1)*(N-S_e)))-(1-C_aditioner)*(-p(8)*S_e.*I_e)).^keta)+...
        sum((abs(I_e_d-C_aditioner*(p(1)*(N-S_e)+p(3)*(t-p(2)).*I_e-p(6))-(1-C_aditioner)*(+p(8)*S_e.*I_e-p(7)*I_e-p(9)*I_e))).^keta)+...
        sum((abs(Q_e_d-C_aditioner*(+p(4)*Q_e+p(5)*I_e+p(6))-(1-C_aditioner)*(+p(7)*I_e-p(10)*Q_e-p(11)*Q_e))).^keta)+...
        sum((abs(R_e_d-C_aditioner*(-p(5)*I_e-(t-p(2))*p(3).*I_e)-(1-C_aditioner)*(+p(9)*I_e+p(10)*Q_e))).^keta)+...
        sum((abs(D_e_d-C_aditioner*(-p(4)*Q_e)-(1-C_aditioner)*(+p(11)*Q_e))).^keta))/N;

A = [];
b = [];
Aeq = [];
beq = [];

l_lub=-30;
h_hub=30;
lb = [l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub,l_lub];
ub = [h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub,h_hub];
options = optimoptions(@fmincon, 'Algorithm' , 'interior-point','Display','off' );
options.Algorithm = 'sqp' ;
% opts.Algorithm = 'dual-simplex' ;
opts.Algorithm = 'interior-point-legacy' ;

[p,fval]=fmincon(obj,p0,A,b,Aeq,beq,lb,ub,[],options);
% [p,fval]= fminsearch(obj,p0);

P_record(iter,:)=p;
F_record(iter)=fval;
end

[M,ind] = min(F_record);
p_estim=P_record(ind,7:11);

parr_record=[parr_record;p_estim];

disp("Error:")
disp(M)

end

data_f=[parr_record pol];
namer="old_siqrd_pol.csv";

writematrix(data_f,namer);


