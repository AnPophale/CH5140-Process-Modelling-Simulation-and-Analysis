%% Case 2 : Fixed F and MR

clc; clear all; close all; format compact; format shortg;

%Global Variables
global xaf ye xe Bh B2 tau Da2
xaf = 0.9;
ye = 0.2;
xe = 0.8;
Bh = 14;
B2 = 4.2;

%% Bifurcation Analysis without Delay
figure;
z_vec = [0.05:0.01:0.9];
Da2_vec = [0.2, 0.4, 0.6, 0.8];
for i = 1:numel(Da2_vec)
    Da2 = Da2_vec(i);
    for j = 1:numel(z_vec)
        f_vec(j) = SS(z_vec(j));
    end
    subplot(2,2,i), plot(z_vec,f_vec,z_vec,0*z_vec)
    xlabel('z')
    ylabel('f(z)')
    title('Da =' + string(Da2))
end
sgtitle('Graphical Method to find number of Steady States');

%From the above plots, we can see that there is only one steady state for
%any value of Da2

%Varying Da2 from 0.2 to 0.7 and finding corresponding SS concentration
%and eigen values
Da2_vec = [0.2:0.001:0.7];
for i = 1:numel(Da2_vec)
    Da2 = Da2_vec(i);
    zSSvec(i) = fsolve(@SS, 0.5, optimoptions('fsolve','Display','none'));
    zSS = zSSvec(i);
    TSS = log(((xaf-zSS) + ((zSS-ye)/(xe-ye))*(xe-xaf))/(Da2*zSS));

    evec(i,1) = max(real(eigfun(zSS,TSS)));
    evec(i,2) = min(real(eigfun(zSS,TSS)));
end

%Plotting steady state concentration vs Da2
figure;
plot(Da2_vec,zSSvec)
xlabel('Da2')
ylabel('z')
title('Steady State concentration vs Da2')

%Plotting eigen values for linear stability analysis vs Da2
figure;
plot(Da2_vec,evec,Da2_vec,0*Da2_vec)
ylim([-2,0.1])
xlabel('Da2');
ylabel('Eigen Value');
title('Eigen Values vs Da2')

%From this plot we can see that the steady state is unstable for Da2 from
%0.433 to 0.511

%% Dynamic analysis without delay
Da2_vec = [0.45, 0.525];
zT0 = [0.2,0.5]; %Initial Conditions
t_final = 150; %Simulation Time

figure;
for i = 1:numel(Da2_vec)
    Da2 = Da2_vec(i);
    nodelay_dyn = ode23s(@dyn, [0,t_final], zT0);
    t = nodelay_dyn.x;
    z = nodelay_dyn.y(1,:);
    subplot(2,1,i), plot(t,z)
    xlabel('t^*');
    ylabel('z');
    title('Da = ' + string(Da2))
end
sgtitle('Concentration vs Time for different Da2 without delay');

%% Phase Plane Plots without delay
zi = [0.1:0.1:0.5]; Ti = [0.5:0.5:3]; %Different initial conditions for the phase plane plot
t_final = 50;

Da2 = 0.45;
figure;
for i = 1:numel(zi)
    for j = 1:numel(Ti)
        nodelay_dyn = ode23s(@dyn, [0,t_final], [zi(i),Ti(j)]);
        z = nodelay_dyn.y(1,:);
        T = nodelay_dyn.y(2,:);
        plot(z,T)
        hold on
        plot(zi(i),Ti(j),'o',...
            'MarkerSize',7)
        hold on
    end
end
xlabel('z');
ylabel('Theta')
title('Phase Plane Plot for Da2=0.45 with no delay');

Da2 = 0.525;
figure;
for i = 1:numel(zi)
    for j = 1:numel(Ti)
        nodelay_dyn = ode23s(@dyn, [0,t_final], [zi(i),Ti(j)]);
        z = nodelay_dyn.y(1,:);
        T = nodelay_dyn.y(2,:);
        plot(z,T)
        hold on
        plot(zi(i),Ti(j),'o',...
            'MarkerSize',7)
        hold on
    end
end
xlabel('z');
ylabel('Theta')
title('Phase Plane Plot for Da2=0.525 with no delay');

%The same result as from the dynamic analysis can also be seen from the
%phase plane plots

%% Bifurcation analysis with delay
%Instead of setting some delay and then finding the corresponding Da2 on the
%critical curve we are varying Da2 and then calcuating tau critical as it is easier to plot

 Da_vec = [0.40:0.001:0.52]; %Varying Da2 to find the corresponding tau critical
 
 for i = 1:numel(Da_vec) 
     Da2 = Da_vec(i);
     sol_critical = fsolve(@delaySS, [0.25,0.05,3],optimoptions('fsolve','Display','None'));
     if sol_critical(2) < 0 %delay should not be negative
         tau_vec(i) = 0;
     else
        tau_vec(i) = sol_critical(2);
     end
 end
figure;
plot(Da_vec,tau_vec)
ylim([0,0.125]);
xlabel('Da2');
ylabel('tau critical');
title('Critical Delay vs Da2 for small values of delay');

%Below the critical curve, the region is unstable, above the curve it is
%stable

%% Dynamic analysis with delay
Da2 = 0.45; %Analysis is done for Da2 = 0.45 and then will be compared to the case without delay
tau_vec = [0.5, 1.5, 2.5, 4]; %Varying the lag to study its effect on the stability

zT0 = [0.2,0.5]; %Initial Conditions
t_final = 150; %Simulation Time

%Calculating the plotting the dynamic response
figure;
for i = 1:numel(tau_vec)
    tau = tau_vec(i);
    delay_dyn = dde23(@delaydyn, tau, zT0,[0,t_final]); %dde23 is used to solve the system of delay differential equations
    t = delay_dyn.x;
    z = delay_dyn.y(1,:);
    subplot(2,2,i), plot(t,z)
    xlabel('t^*')
    ylabel('z')
    title('Tau = ' + string(tau))
    hold on
end

sgtitle('Concentration vs Time for different Da2 with delay')

%For no delay, Da2 = 0.45 always gives an unstable steady state.
%Depending on the delay, we see that there can be a stable steady
%state for the same Da1 but with delay.

%% Function Definitions
%Function to solve for steady state concentration without delay
function g = SS(z)
    global xaf ye xe Bh B2 Da2
    T = log(((xaf-z) + ((z-ye)/(xe-ye))*(xe-xaf))/(Da2*z));
    g = -T*(1+B2) + Bh*Da2*z*exp(T);
end

%Function to find Jacobain of the system of equations without delay
function f = eigfun(z,T)
    global xaf ye xe Bh B2 Da2
    J = [-1 + ((xe-xaf)/(xe-ye)) - Da2*exp(T), -Da2*z*exp(T);
        Bh*Da2*exp(T), -(1+B2) + Bh*Da2*z*exp(T)];
    f = eig(J);
end

%Function for dynamic analysis without delay
function f = dyn(t,zT)
    global xaf ye xe Bh B2 Da2
    
    z = zT(1);
    T = zT(2);
    
    f(1,1) = (((xe-z)*(xaf-ye))/(xe-ye)) - Da2*z*exp(T);
    f(2,1) = -T*(1+B2) + Bh*Da2*z*exp(T);
end

%Function to find the tau critical for a given Da2 for the system with delay
function f = delaySS(X)
 global xaf ye xe Bh B2 Da2
    z = X(1);
    tau = X(2);
    w = X(3);

    T = log(((xe-z)*(xaf-ye))/((xe-ye)*Da2*z));
    f(1) = -T*(1+B2) + Bh*Da2*z*exp(T);
    
    f(2) = -(w^2) -w*(((xe-xaf)/(xe-ye))*sin(w*tau))+ ...
        (1 + Da2*exp(T) - ((xe-xaf)/(xe-ye))*cos(w*tau))*(1+B2 - Bh*Da2*z*exp(T))+ ...
        Bh*(Da2^2)*z*exp(2*T);
    
    f(3) = w*(2+B2 + Da2*exp(T) - Bh*Da2*z*exp(T) - ((xe-xaf)/(xe-ye))*cos(w*tau))+ ...
        ((xe-xaf)/(xe-ye))*sin(w*tau)*(1+B2 - Bh*Da2*z*exp(T));
end

%Function to solve the system of delay differential equations
%For dynamic analysis with delay
function f = delaydyn(t,zT,y)
 global xaf ye xe Bh B2 Da2
    ylag = y(1);

    z = zT(1);
    T = zT(2);

    f(1,1) = (xaf-z) + ((ylag - ye)/(xe - ye))*(xe - xaf) - Da2*z*exp(T);
    f(2,1) = -T*(1+B2)+ Bh*Da2*z*exp(T);  
end
