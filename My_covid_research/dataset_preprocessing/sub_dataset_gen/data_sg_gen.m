clear all
close all
clc

len_data=573;

num_par=9;

for i=0:len_data
    tmp_name="D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\input\states_";
    curr_stat=csvread(join([tmp_name,num2str(i),".csv"]));
    name_input=join(["D:\Research_work\covid_parameters_2\siqrd\sub_dataset_gen\input\states_",num2str(i),".csv"]);
    writematrix(curr_stat,name_input);
end

t=0:num_par-1;
fs=3;
t_new=linspace(0,num_par-1,fs*length(t));

for i=0:len_data
    tmp_name="D:\Research_work\My_covid_research\dataset_preprocessing\parameter_estim\output\parameters_";
    curr_para=csvread(join([tmp_name,num2str(i),".csv"]));
    p=curr_para(1:num_par);
    p_new = spline(t,p,t_new);
    name_output=join(["output_norm_interpolated\parameters_",num2str(i),".csv"]);
    writematrix(p_new,name_output);
end

for i=0:len_data
    tmp_name="output_norm_interpolated\parameters_";
    curr_para2=csvread(join([tmp_name,num2str(i),".csv"]));
    order = 1;
    framelen = 3;
    curr_para2=curr_para2-mean(curr_para2);
    p_sg = sgolayfilt(curr_para2+1*rand(1,length(t_new)),order,framelen);
%     p_sg = sgolayfilt(curr_para2+0*rand(1,length(t_new)),order,framelen);
    p_para= polyfit(t_new,p_sg,5);
    p_sg=polyval(p_para,t_new);
    p_sg=tanh(0.05*p_sg);
    name_output2=join(["output_sg_interpolated\parameters_",num2str(i),".csv"]);
    writematrix(p_sg,name_output2);
end


figure
subplot(3,1,1)
plot(t,p)
subplot(3,1,2)
plot(t_new,p_new)
subplot(3,1,3)
plot(t_new,p_sg)
