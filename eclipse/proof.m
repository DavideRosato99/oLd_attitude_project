close all; clc; 

a = [-30:0.01:30] *pi/180;
D = 200000+6371;
Wf = [];
Ys = [];


for alpha = a
    [Wf1] = SRP_eclipse_model(alpha, D);
    Wf = [Wf Wf1];
    Ys1 = D*tan(alpha);
    Ys = [Ys Ys1];
end


plot(a,Wf)


