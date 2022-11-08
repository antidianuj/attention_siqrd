clear all
close all
clc

my_ind=0;
observe=3;

err_ind=0;
offsetter=173;

hor=173;
curr_par=csvread("D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\data_augmentation\parameter_generation\parameter_generated.csv");
a=size(curr_par);
% Adjust C between 0 and 1

threshold=200;



for i=1:a(1)
    for j=1:hor
        if my_ind<threshold
        curr_stat=csvread(join(["D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\input\states_",num2str(j-1),".csv"]));

        days=linspace(0,observe-1,observe);
        Fs=20; 
        t=linspace(0,length(days)-1,length(days)*Fs);
        my_time=t;
        
        initcond=curr_stat(1,2:6);
        N=sum(initcond);




        zeta=curr_par(i,1);
        beta1=curr_par(i,2);
        beta2=curr_par(i,3);
        alpha1=curr_par(i,4);
        alpha2=curr_par(i,5);
        delta=curr_par(i,6);
        K=curr_par(i,7);
        to=curr_par(i,8);
        gamma=curr_par(i,9);
        
        
        SIQRD = @(my_time,x) ([ -zeta*(N-x(1))+gamma*x(4)
                zeta*(N-x(1))-beta1*(to-my_time)*x(2)-beta2*x(2)-K-delta*x(2)
                beta2*x(2)+K-alpha2*x(3)-alpha1*x(3)
                beta1*(to-my_time)*x(2)+alpha1*x(3)-gamma*x(4)
                delta*x(2)*x(5)+alpha2*log(1+x(3))]);



        % opts = odeset( 'RelTol' ,1e-1, 'maxstep' ,1e-2);
        
        [my_time,y] = ode15s(SIQRD, my_time,initcond);
        
        
        
        if any(any(isnan(abs(y))==1)) || any(any(isinf(abs(y))==1))  || any(any(abs(y)-N>20))
            err_ind=err_ind+1;
            disp(join(["-----------",num2str(err_ind),"-------------"]));
            disp(i-1);
            disp("---fuck----");
            continue
        
        else
        S=y(1:Fs:end,1);
        I=y(1:Fs:end,2);
        Q=y(1:Fs:end,3);
        R=y(1:Fs:end,4);
        D=y(1:Fs:end,5);
        
        
        
        try
        yprime=[days' S I Q R D];
%         plot(days',S)
%         hold on
%         plot(days',I)
        para_new=[zeta,beta1,beta2,alpha1,alpha2,delta,K,to,gamma];
        name_output=join(["D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\output\parameters_",num2str(offsetter+1+my_ind),".csv"]);
        name_input=join(["D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\input\states_",num2str(offsetter+1+my_ind),".csv"]);
        writematrix(yprime,name_input);
        writematrix(para_new,name_output);
        disp("----------------")
        disp(offsetter+1+my_ind)
        my_ind=my_ind+1;


        catch ME
            disp("Oops")
            continue


        end
end
    end



end
end
