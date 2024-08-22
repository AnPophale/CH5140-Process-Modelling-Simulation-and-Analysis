%% Case 1 : Fixed F0 and MR

clc; clear all; close all; format compact; format shortg;

%Global Variables
global xaf ye xe Bh B1 tau Da1
xaf = 1;
ye = 0.1;
xe = 0.9;
Bh = 8.8;
B1 = 4;

%% Bifurcation Analysis without Delay
figure;
z_vec = [0.05:0.01:0.9];
Da1_vec = [0.4, 0.457, 0.571, 0.6];
for i = 1:numel(Da1_vec)
    Da1 = Da1_vec(i);
    for j = 1:numel(z_vec)
        f_vec(j) = SS(z_vec(j));
    end
    subplot(2,2,i), plot(z_vec,f_vec,z_vec,0*z_vec)
    xlabel('z')
    ylabel('f(z)')
    title('Da =' + string(Da1))
    ylim([-50,10]);
end
sgtitle('Graphical Method to find number of Steady States');

%From the above plots, we can see that for Da1 < 0.457 there is no steady
%state, at 0.457 we have on SS and beyond it there are two steady states

%Varying Da1 from 0.457 to 0.7 and finding corresponding SS concentration
%and eigen values
Da1_vec = [0.457:0.001:0.7];
for j = 1:numel(Da1_vec)
    Da1 = Da1_vec(j);
    zSS1vec(j) = fsolve(@SS,0.4,optimoptions('fsolve','Display','None'));
    zSS2vec(j) = fsolve(@SS,0.8,optimoptions('fsolve','Display','None'));

    zSS1 = zSS1vec(j); zSS2 = zSS2vec(j);
    TSS1 = log((xaf - ye)/(zSS1*Da1));
    TSS2 = log((xaf - ye)/(zSS2*Da1));

    evec1(j,1) = max(real(eigfun(zSS1, TSS1)));
    evec1(j,2) = min(real(eigfun(zSS1, TSS1)));

    evec2(j,1) = max(eigfun(zSS2, TSS2));
    evec2(j,2) = min(eigfun(zSS2, TSS2));
end    

%Plotting steady state concentration vs Da1
figure;
plot(Da1_vec, zSS1vec, '-r', Da1_vec, zSS2vec, '-b')
title('Steady State concentration vs Da1');
xlabel('Da1');
ylabel('z');
legend('SS1','SS2')

%Plotting eigen values for linear stability analysis vs Da1
figure;
subplot(1,2,1), plot(Da1_vec,evec1,'-o')
hold on
xlabel('Da1');
ylabel('Eigen Value');
title('SS1 (Lower)')

subplot(1,2,2), plot(Da1_vec,evec2,'-o')
xlabel('Da1');
ylabel('Eigen Value');
title('SS2 (Higher)')

sgtitle('Eigen values corresponding to the SS vs Da1');

%From these plots we can see that the higher SS is always unstable while
%the lower is unstable only when Da1 is between 0.464 and 0.571.

%Hence there are no stable steady states for Da1 lying between 0.464 and
%0.571. For other values of Da1, only the lower SS is stable.

%% Dynamic Analysis without delay
Da1_vec = [0.5, 0.6]; %Analysis is done for Da1 = 0.5, 0.6 to see the change in stability of SS
zT0 = [0.2,0.5]; %Initial Conditions
t_final = 100; %Simulation Time

figure;
for i = 1:numel(Da1_vec)
    Da1 = Da1_vec(i);
    nodelay_dyn = ode23s(@dyn, [0,t_final], zT0);
    t = nodelay_dyn.x;
    z = nodelay_dyn.y(1,:);
    subplot(2,1,i), plot(t,z)
    xlabel('t^*');
    ylabel('z');
    title('Da = ' + string(Da1))
end

sgtitle('Concentration vs Time for different Da1 with no delay');

%As expected for Da1 = 0.5 we get sustained oscillations but for Da1 = 0.6,
%the system reaches the lower steady state.

%% Phase Plane Plots without delay
zi = [0.1:0.1:0.5]; Ti = [0.5:0.5:3]; %Different initial conditions for the phase plane plot
t_final = 50;

Da1 = 0.5;
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
title('Phase Plane Plot for Da1=0.5 with no delay');

Da1 = 0.6;
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
title('Phase Plane Plot for Da1=0.6 with no delay');

%The same result as from the dynamic analysis can also be seen from the
%phase plane plots

%% Bifurcation analysis with delay
%Instead of setting some delay and then finding the corresponding Da1 on the
%critical curve we are varying Da1 and then calcuating tau critical as it is easier to plot

 Da_vec = [0.45:0.001:0.58]; %Varying Da1 to find the corresponding tau critical
 
 for i = 1:numel(Da_vec) 
     Da1 = Da_vec(i);
     sol_critical = fsolve(@delaySS, [0.5,0.5,2],optimoptions('fsolve','Display','None'));
     if sol_critical(2) < 0 %delay should not be negative
         tau_vec(i) = 0;
     else
        tau_vec(i) = sol_critical(2);
     end
 end
figure;
plot(Da_vec,tau_vec)
ylim([0,0.125]);
xlabel('Da1');
ylabel('tau critical');
title('Critical Delay vs Da1 for small values of delay');

%Below the critical curve, the region is unstable, above the curve it is
%stable

%% Dynamic analysis with delay
Da1 = 0.5; %Analysis is done for Da1 = 0.5 and then will be compared to the case without delay
tau_vec = [0.05, 0.125, 2.5, 3.75, 5, 8]; %Varying the lag to study its effect on the stability

zT0 = [0.2,0.5]; %Initial Conditions
t_final = 200; %Simulation Time

%Calculating the plotting the dynamic response
figure;
for i = 1:numel(tau_vec)
    tau = tau_vec(i);
    delay_dyn = dde23(@delaydyn, tau, zT0,[0,t_final]); %dde23 is used to solve the system of delay differential equations
    t = delay_dyn.x;
    z = delay_dyn.y(1,:);
    subplot(3,2,i), plot(t,z)
    xlabel('t^*')
    ylabel('z')
    title('Tau = ' + string(tau))
    hold on
end

sgtitle('Concentration vs Time for different Da1 with delay')

%For no delay, Da1 = 0.5 always gives an unstable steady state.
%Depending on the delay, we see that there can be a stable steady
%state for the same Da1 but with delay.

%% Function Definitions
%Function to solve for steady state concentration without delay
function g = SS(z)
    global xaf ye xe Bh B1 Da1
    g = -log((xaf - ye)/(Da1*z))*(1 + B1+ (z - ye)/(xe - z)) +  Bh*Da1*z*((xaf - ye)/(Da1*z));
end

%Function to find Jacobain of the system of equations without delay
function f = eigfun(z,T)
    global ye xe Bh B1 Da1
    J = [-Da1*exp(T), -Da1*z*exp(T);
    T*(ye-xe)/((xe-z)^2) + Bh*Da1*exp(T),  -(1+B1 +((z-ye)/(xe-z))) + Bh*Da1*z*exp(T)];
    f = eig(J);
end

%Function for dynamic analysis without delay
function f = dyn(t,zT)
    global xaf ye xe Bh B1 Da1
    
    z = zT(1);
    T = zT(2);
    
    f(1,1) = xaf - ye - Da1*z*exp(T);
    f(2,1) = -T*(1+B1 + ((z-ye)/(xe-z))) + Bh*Da1*z*exp(T);
end

%Function to find the tau critical for a given Da1 for the system with delay
function f = delaySS(X)
    global xaf ye xe Bh B1 Da1 
    z = X(1);
    tau = X(2);
    w = X(3);

    T = log((xaf-ye)/(Da1*z));
    f(1) = -T*(1+B1+((z-ye)/(xe-z))) + Bh*Da1*z*exp(T);
    
    f(2) = -(w^2) + w*((ye-xe)/(xe-z))*sin(w*tau)+ ...
        (((ye-xe)/(xe-z))*(cos(w*tau) -1) +Da1*exp(T))*(B1 + ((xe-ye)/(xe-z)) - Bh*Da1*z*exp(T))- ...
        Da1*exp(T)*z*((T*(xe-ye)/((xe-z)^2))*cos(w*tau) - Bh*Da1*exp(T));
    
    f(3) = w*(  B1 + ((xe-ye)/(xe-z)) - Bh*Da1*z*exp(T) + ((ye-xe)/(xe-z))*(cos(w*tau) -1) + Da1*exp(T))- ...
        ((ye-xe)/(xe-z))*sin(w*tau)*(B1 + ((xe-ye)/(xe-z)) -Bh*Da1*z*exp(T))+ ...
        Da1*exp(T)*z*T*((xe-ye)/((xe-z)^2))*sin(w*tau);
end

%Function to solve the system of delay differential equations
%For dynamic analysis with delay
function f = delaydyn(t,zT,y)
 global xaf ye xe Bh B1 Da1
    ylag = y(1);

    z = zT(1);
    T = zT(2);

    f(1,1) = (xaf-z) + ((ylag - ye)/(xe - ylag))*(xe - z) - Da1*z*exp(T);
    f(2,1) = -T*(1+B1+ (ylag - ye)/(xe - ylag)) + Bh*Da1*z*exp(T);  
end
