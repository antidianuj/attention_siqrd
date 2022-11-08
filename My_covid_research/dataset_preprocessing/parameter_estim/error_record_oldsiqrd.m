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
pol=data(3:342,6:6+10-1);


par_len=5;

N=mean(S+I+Q+R+D);

l_h=50+10+150;
h_h=300;

Sp=S(l_h:h_h);
Ip=I(l_h:h_h);
Qp=Q(l_h:h_h);
Rp=R(l_h:h_h);
Dp=D(l_h:h_h);


observ_len=3;

for observe=1:h_h-l_h+1-observ_len+1


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
tx=t;
t=t(1:end-1);


lamb1=0.0;
lamb2=0.0000;



tries=70;
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
        ];
%     p0 = (200)*rand(par_len,1);
keta=2;
E=1;
C_aditioner=0;

% alpha=p1
% beta=p2
% gamma1=p3
% gamma2=p4
% delta=p5

obj=@(p)(sum((abs(S_e_d-(1-C_aditioner)*(-p(2)*S_e.*I_e))).^keta)+...
        sum((abs(I_e_d-(1-C_aditioner)*(+p(2)*S_e.*I_e-p(1)*I_e-p(3)*I_e))).^keta)+...
        sum((abs(Q_e_d-(1-C_aditioner)*(+p(1)*I_e-p(4)*Q_e-p(5)*Q_e))).^keta)+...
        sum((abs(R_e_d-(1-C_aditioner)*(+p(3)*I_e+p(4)*Q_e))).^keta)+...
        sum((abs(D_e_d-(1-C_aditioner)*(+p(5)*Q_e))).^keta))/N;

A = [];
b = [];
Aeq = [];
beq = [];

l_lub=-30;
h_hub=30;
lb = [l_lub,l_lub,l_lub,l_lub,l_lub];
ub = [h_hub,h_hub,h_hub,h_hub,h_hub];
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

p=P_record(ind,:);

alpha=p(1);
beta=p(2);
gamma1=p(3);
gamma2=p(4);
delta=p(5);


nSIQRD = @(tx,x) ([ -beta*x(1)*x(2)
        beta*x(1)*x(2)-alpha*x(2)-gamma1*x(3)
        alpha*x(2)-gamma2*x(3)-delta*x(3)
        gamma1*x(2)+gamma2*x(3)
        delta*x(3)]);


initcond=[S_e(1) I_e(1) Q_e(1) R_e(1) D_e(1)];
% opts = odeset( 'RelTol' ,1e-1, 'maxstep' ,1e-2);       
[tx,y] = ode15s(nSIQRD, tx ,initcond);
        

S_p=y(1:Fs:end,1);
I_p=y(1:Fs:end,2);
Q_p=y(1:Fs:end,3);
R_p=y(1:Fs:end,4);
D_p=y(1:Fs:end,5);


my_error=(sum(abs(S-S_p))+sum(abs(I-I_p))+sum(abs(Q-Q_p))+sum(abs(R-R_p))+sum(abs(D-D_p)))/(5*length(D));
% Mean absolute error
err_record=[err_record;my_error];


end


figure
plot(err_record);
xlabel('Days');
ylabel('MAE');
grid on

figure
numBins = 400;
histogram(err_record, numBins, 'EdgeColor', 'none');
grid on;
xlabel('MAE', 'FontSize', 15);
ylabel('Density', 'FontSize', 15);
% Compute mean and standard deviation.
mu = mean(err_record);
sigma = std(err_record);
% Indicate those on the plot.
% xline(mu, 'Color', 'g', 'LineWidth', 2);
% xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
% xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
% ylim([0, 20]); % Give some headroom above the bars.
% xlim([0, 40]); % Give some headroom above the bars.
yl = ylim;
sMean = sprintf('  Mean = %.9f\n  SD = %.9f', mu, sigma);
% Position the text 90% of the way from bottom to top.
text(mu/2, 0.9*yl(2), sMean, 'Color', 'r', ...
	'FontWeight', 'bold', 'FontSize', 12, ...
	'EdgeColor', 'b');
sMean2= sprintf('Mean = %.3f.  SD = %.3f', numBins, mu, sigma);
% title(sMean2, 'FontSize', 15);


