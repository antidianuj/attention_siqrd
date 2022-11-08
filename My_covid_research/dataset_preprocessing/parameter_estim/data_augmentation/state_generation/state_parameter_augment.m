
clear all
close all
clc

thresher=8500;
offseter=373;
Lamb=1;

threshold=200;

data=csvread("state_generated.csv");
err_record=[];
parr_record=[];
observer=0;


S=data(:,1)*Lamb;
I=data(:,2)*Lamb;
Q=data(:,3)*Lamb;
R=data(:,4)*Lamb;
D=data(:,5)*Lamb;

par_len=9;



l_h=1;
h_h=length(data);


Sp=S(l_h:h_h);
Ip=I(l_h:h_h);
Qp=Q(l_h:h_h);
Rp=R(l_h:h_h);
Dp=D(l_h:h_h);


observ_len=3;

for observe=1:ceil((h_h-l_h)/observ_len)-1

if observer<threshold
S=Sp((observe-1)*observ_len+1:observe*observ_len);
I=Ip((observe-1)*observ_len+1:observe*observ_len);
Q=Qp((observe-1)*observ_len+1:observe*observ_len);
R=Rp((observe-1)*observ_len+1:observe*observ_len);
D=Dp((observe-1)*observ_len+1:observe*observ_len);
days=linspace(0,observ_len-1,observ_len);

order = 1;
framelen = 11;

plot(S)
hold on
plot(I)
plot(Q)
plot(R)
plot(D)
grid on

if min(S)<0
    S=S+abs(min(S));
end

if min(I)<0
    I=I+abs(min(I));
end


if min(Q)<0
    Q=Q+abs(min(Q));
end


if min(R)<0
    R=R+abs(min(R));
end


if min(D)<0
    D=D+abs(min(D));
end






Fs=20; 
t=linspace(0,length(days)-1,length(days)*Fs);

S_e=spline(days,S,t);
I_e=spline(days,I,t);
Q_e=spline(days,Q,t);
R_e=spline(days,R,t);
D_e=spline(days,D,t);

if min(S_e)<0
    S_e= sgolayfilt(S_e,order,framelen);
    S_e=S_e+abs(min(S_e));
    
end

if min(I_e)<0
    I_e= sgolayfilt(I_e,order,framelen);
    I_e=I_e+abs(min(I_e));
    
end


if min(Q_e)<0
    Q_e= sgolayfilt(Q_e,order,framelen);
    Q_e=Q_e+abs(min(Q_e));
    
end


if min(R_e)<0
    R_e= sgolayfilt(R_e,order,framelen);
    R_e=R_e+abs(min(R_e));
end


if min(D_e)<0
    D_e= sgolayfilt(D_e,order,framelen);
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
% [p,fval]= fminsearch(obj,p0);

P_record(iter,:)=p;
F_record(iter)=fval;
end

[M,ind] = min(F_record);
p_estim=P_record(ind,:);

err_record=[err_record, M];
parr_record=[parr_record;p_estim];
disp(join([num2str(observe),": ",num2str(M)]))

% if(M<5)
if(M<thresher)
try
    
    data_in=[days',S,I,Q,R,D];
   
    data_out=p_estim;
    namer_in=join(["D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\input\states_",num2str(observer+offseter+1),".csv"]);
    namer_out=join(["D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\output\parameters_",num2str(observer+offseter+1),".csv"]);
    writematrix(data_in,namer_in);
    writematrix(data_out,namer_out);
    observer=observer+1;
    disp("************")
    disp(observer+offseter+1)
    
    
catch ME
    continue

end
end

end
end
title("Example Suspected and Infected")

figure
plot(err_record);
grid on
title("Error Record")

% figure
% plot(parr_record(:,:));
% grid on
% title("Parameter Record")
% legend("alpha1","t1","gamma11","gamma21","beta1","K","alpha2","beta2","gamma12","gamma22","delta2","C");
% 
% disp(max(err_record))
% 
% curr_par=p_estim;
% min=0;
% maxx=0;
% alpha=curr_par(1)+min+(maxx-min)*rand;
% t1=curr_par(2)+min+(maxx-min)*rand;
% gamma1=curr_par(3)+min+(maxx-min)*rand;
% gamma2= curr_par(4)+min+(maxx-min)*rand;
% beta= curr_par(5)+min+(maxx-min)*rand;
% K= curr_par(6)+min+(maxx-min)*rand;
% alpha2=curr_par(7)+min+(maxx-min)*rand;
% beta2=curr_par(8)+min+(maxx-min)*rand;
% gamma12=curr_par(9)+min+(maxx-min)*rand;
% gamma22=curr_par(10)+min+(maxx-min)*rand;
% delta2=curr_par(11)+min+(maxx-min)*rand;
% C=curr_par(12)+min+(maxx-min)*rand;
%  
% days=linspace(0,observ_len-1,observ_len);
% Fs=20; 
% t=linspace(0,length(days)-1,length(days)*Fs);
% my_time=t;
% 
% 
% SIQRD = @(my_time,x) ([ E*C*(alpha*(N-x(1)))+(1-C)*(-beta2*x(1)*x(2))
%             E*C*(-alpha*(N-x(1))-(my_time-t1)*gamma1*x(2)-K)+(1-C)*(+beta2*x(1)*x(2)-alpha2*x(2)-gamma12*x(2))
%             E*C*(-gamma2*x(3)-beta*x(2)+K)+(1-C)*(+alpha2*x(2)-gamma22*x(3)-delta2*x(3))
%             E*C*(+(my_time-t1)*gamma1*x(2)+beta*x(2))+(1-C)*(+gamma12*x(2)+gamma22*x(3))
%             E*C*(gamma2*x(3))+(1-C)*(+delta2*x(3))]);
% 
% 
% %     opts = odeset( 'RelTol' ,1e-1, 'maxstep' ,1e-2);
% initcond=data_in(1,2:6);
% [my_time,y] = ode15s(SIQRD, my_time, initcond);
% 
% S=y(1:Fs:end,1);
% I=y(1:Fs:end,2);
% Q=y(1:Fs:end,3);
% R=y(1:Fs:end,4);
% D=y(1:Fs:end,5);
% 
% figure
% plot(days,I,days,Q,days,R,days,D);
% figure
% plot(data_in(:,1),data_in(:,3),data_in(:,1),data_in(:,4),data_in(:,1),data_in(:,5),data_in(:,1),data_in(:,6));