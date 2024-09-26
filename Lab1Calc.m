clear
close all
clc
%line 37 and 64 to stop 
% Adding Starting var.
h = struct2cell(load('H.mat'));
w = h{1};
px2m = (1/37)*(0.0254);
Pb_H = (1024-w{1}).*px2m;
PPb_H = (1024-w{2}).*px2m;
Sb_H = (1024-w{3}).*px2m;
Tb_H = (1024-w{4}).*px2m;
t1 = [1:length(Pb_H)]./200;
t2 = [1:length(PPb_H)]./200;
t3 = [1:length(Sb_H)]./200;
t4 = [1:length(Tb_H)]./200;
Pb_mass = 26/1000;
PPb_mass = 3/1000;
Sb_mass = 27/1000;
Tb_mass = 57/1000;

% Q1 PLotting DvT Graphs for al four curves
hold on
plot(t1 , Pb_H,'b-', 'linewidth', 2);
plot(t2 , PPb_H,'r-', 'linewidth', 2);
plot(t3 , Sb_H,'y-', 'linewidth', 2);
plot(t4 , Tb_H,'g-', 'linewidth', 2);
hold off
title('Height of balls per Time')
xlabel('Time (s)') 
ylabel('Hieght (m)')
legend({'Pickle Ball','Ping Pong Ball','Squash Ball','Tennis Ball'},'Location','Northeast')
    shg
hold off

%Q2 Find the and plot the VvT graph of Pb_H at three different time step interval
V = 0;
tp = t1(1);
res = 1;
for k = 2:length(Pb_H)/res
    V(k) = (Pb_H((k*res)-1)-Pb_H(k*res))/(t1((res*k)-1)-t1(k*res));
    tp(k) = t1(k*res);
end    

plot(tp , V,'b-', 'linewidth', 2);
hold on
V = 0;
tp = t1(1);
res = 10;
for k = 2:length(Pb_H)/res
    V(k) = (Pb_H((k*res)-1)-Pb_H(k*res))/(t1((res*k)-1)-t1(k*res));
    tp(k) = t1(k*res);
end    

plot(tp , V,'r-', 'linewidth', 2);
V = 0;
tp = t1(1);
res = 25;
for k = 2:length(Pb_H)/res
    V(k) = (Pb_H((k*res)-1)-Pb_H(k*res))/(t1((res*k)-1)-t1(k*res));
    tp(k) = t1(k*res);
end    
plot(tp , V,'g-', 'linewidth', 2);
title('Velocity Vs Time at 3 Temporal Spacings')
xlabel('Time (s)') 
ylabel('Velocity (m/s)')
legend({'Resolution 1','Resolution 10','Resolution 25'},'Location','Northeast')
hold off

%Q3 Calc CoR in terms of height and velocity for each ball and calc the
%uncertainitees
Nam = {Pb_H ,PPb_H ,Sb_H ,Tb_H};
unc_H = 0.003175
for k = 1:length(Nam)
    H1_0_Max = max(Nam{k});
    [o , low_B] = min(Nam{k});
    H1_1_Max = max(Nam{k}(low_B:end));
    CoR_H(k) = sqrt(H1_1_Max/H1_0_Max);
    Unc_CoR(k) = unc_H*sqrt((2 * H1_0_Max^2 + H1_1_Max^2)/(4*H1_0_Max^3*H1_1_Max))

end

Nam = {Pb_H ,PPb_H ,Sb_H ,Tb_H};
ts = {t1, t2, t3, t4};
for v = 1:length(Nam)
    V = 0;
    tp = ts{v}(1);
    res = 3;
    for k = 2:length(Nam{v})/res
        V(k) = (Nam{v}((k*res)-1)-Nam{v}(k*res))/(ts{v}((res*k)-1)-ts{v}(k*res));
        tp(k) = ts{v}(k*res);
    end    
    [V_max , I1] = max(V);
    [V_min , I2] = min(V);
    T1 = ts{v}(I1*res);
    T2 = ts{v}(I2*res);
    test = (2*(V_min^2-V_max^2))/(T1-T2)
    test1 = (1)/(T1-T2)
    Unc_CoRv(v) = unc_H * V_min^2 * sqrt((2*(V_min^2-V_max^2))/(T1-T2));
    CoR_V(v) = V_max/abs(V_min);

end
CoR_H
Unc_CoR
CoR_V
Unc_CoRv


%Q4 find the chanve in V at point of impact and scale with Mass for each
%ball as force(F), Then plot FvCoR as four points include error bars

b_Mass = {Pb_mass , PPb_mass , Sb_mass , Tb_mass};
for v = 1:length(Nam)
    [ o , Imp ] = min(Nam{v});
    ts = {t1, t2, t3, t4};
    V = 0;
    tp = ts{v}(1);
    res = 1;
    for k = 2:length(Nam{v})/res
        V(k) = (Nam{v}((k*res)-1)-Nam{v}(k*res))/(ts{v}((res*k)-1)-ts{v}(k*res));
        tp(k) = ts{v}(k*res);
    end
    [M1 , I1] = max(V(Imp-5:Imp+5));
    [M2 , I2] = min(V(Imp-5:Imp+5));
    diff_V = M1 - M2;
    diff_T = ts{v}(I1) - ts{v}(I2);
    F = b_Mass{v} * (diff_V/diff_T);
    Unc_F(v) = unc_H * 2 * sqrt(1/(diff_T)^3);
    scatter(F , CoR_H(v) , "Filled")
    hold on
    errorbar( F , CoR_H(v) , Unc_CoR(v));
    errorbar( F , CoR_H(v) , Unc_F(v), 'horizontal');
end
Unc_F
title('Force during Impact vs. Coefficient of Restitution')
axis("padded")
xlabel('Force (N)') 
ylabel('Coeff of Resitution')
legend({'Pickle Ball','','','Ping Pong Ball','','','Squash Ball','','','Tennis Ball','',''},'Location','Northeast')
    shg
