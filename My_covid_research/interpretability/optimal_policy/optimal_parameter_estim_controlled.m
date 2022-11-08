clear all
close all
clc


data=csvread("D:\Research_work\My_covid_research\approximation_model\future_data_est.csv");
data_par=csvread("D:\Research_work\My_covid_research\approximation_model\parameters_est.csv");


init_s_data=csvread("D:\Research_work\COVID_paper\My_covid_research\dataset_preprocessing\parameter_estim\input_new\states_ 0 .csv");

my_init_cond=init_s_data(1,2:end);



zeta=data_par(:,1);
beta1=data_par(:,2);
beta2=data_par(:,3);
alpha1=data_par(:,4);
alpha2=data_par(:,5);
delta=data_par(:,6);
K=data_par(:,7);
to=data_par(:,8);
gamma=data_par(:,9);


I=data(:,2);
D=data(:,5);
S=data(:,1);
Q=data(:,3);
t=linspace(0,15,15)';
N=51343545;
tries=50;
P_record=zeros(tries,2);
F_record=zeros(tries,1);
minner=-10;
maxxer=10;




k=10;
u=0.1;
for iter=1:tries
obj=@(p)(sum(((zeta+p(1)).*(N-S)-beta1.*(to-t).*I-(beta2+p(2)).*I-K-(delta).*I+u*((delta).*I.*D+(alpha2).*log(1+Q))+k*(I+u*D)).^2))/N;
    p0=[
        -minner+maxxer*rand
       -minner+maxxer*rand
        ];


A = [];
b = [];
Aeq = [];
beq = [];

l_lub=-20;
h_hub=20;
lb = [l_lub,l_lub];
ub = [h_hub,h_hub];
options = optimoptions(@fmincon, 'Algorithm' , 'interior-point','Display','off','PlotFcns',{@optimplotfval,});
options.Algorithm = 'sqp' ;
% opts.Algorithm = 'dual-simplex' ;
opts.Algorithm = 'interior-point-legacy' ;

[p,fval]=fmincon(obj,p0,A,b,Aeq,beq,lb,ub,[],options);
P_record(iter,:)=p;
F_record(iter)=fval;
end

[M,ind] = min(F_record);
peta=P_record(ind,:);



% zeta=data_par(:,1);
% beta1=data_par(:,2);
% beta2=data_par(:,3);
% alpha1=data_par(:,4);
% alpha2=data_par(:,5);
% delta=data_par(:,6);
% K=data_par(:,7);
% to=data_par(:,8);
% gamma=data_par(:,9);




zeta=mean(data_par(:,1))+peta(1);
beta1=mean(data_par(:,2));
beta2=mean(data_par(:,3))+peta(2);
alpha1=mean(data_par(:,4));
alpha2=mean(data_par(:,5));
delta=mean(data_par(:,6));
K=mean(data_par(:,7));
to=mean(data_par(:,8));
gamma=mean(data_par(:,9));

Fs=20;
my_time=linspace(0,15-1,15*Fs);
SIQRD = @(my_time,x) ([ -zeta*(N-x(1))+gamma*x(4)
                zeta*(N-x(1))-beta1*(to-my_time)*x(2)-beta2*x(2)-K-delta*x(2)
                beta2*x(2)+K-alpha2*x(3)-alpha1*x(3)
                beta1*(to-my_time)*x(2)+alpha1*x(3)-gamma*x(4)
                delta*x(2)*x(5)+alpha2*log(1+x(3))]);



opts = odeset( 'RelTol' ,1e-1, 'maxstep' ,1e-2);
        
[my_time,y] = ode15s(SIQRD, my_time,my_init_cond);


S=y(1:Fs:end,1);
I=y(1:Fs:end,2);
Q=y(1:Fs:end,3);
R=y(1:Fs:end,4);
D=y(1:Fs:end,5);

o_states=[S,I,Q,R,D];

writematrix(o_states,'optimal_state_trajecotry.csv');


par=[zeta,beta1,beta2,alpha1,alpha2,delta,K,to,gamma];
namer="optimal_par.csv";

writematrix(par,namer);

par=[mean(data_par(:,1)),beta1,mean(data_par(:,3)),alpha1,alpha2,delta,K,to,gamma];
namer="passive_par.csv";

writematrix(par,namer);

% f = @(x,y)(sum((zeta.*(N-S)+beta1.*(t-to).*I-(K+y)+u*(alpha+x).*Q+k*(I+u*D)).^2))/1e18;
% fsurf(f,[-30,30],'ShowContours','on')




