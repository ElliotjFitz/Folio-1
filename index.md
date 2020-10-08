## Welcome to Elliot's Digital Folio

### Simulation of 5098 sand particles in 2D under axial loading

In this project, I used R and MATLAB to simulate and analyse a system of 5098 sand particles under axial loading in 2D. Here I used a measure of betweenness centrality to predict the shear band's location prior to peak loading stress. The operating assumtion here is that force travels along the path of least resistance. Notice that the colouring of the nodes show the orange nodes as having the highest betweenness centrality measure. This wonderfully coincides with the approximate location of the V shaped shear band.

![Shear Band Sim](https://github.com/ElliotjFitz/Page/blob/gh-pages/0.02CVElliot.avi?raw=true)



### Kalman Filter Implementation in MATLAB
The following is the implementation of a Kalman Filter coded without using ss and kalmf. The scenario presented is as follows: There are 3 dams with only two having working sensors and so we must use a Kalman filter to estimate the missing sensor's value. The timescale here is to the order of hours and the model is linearised. 
```
# Filter code

                                                                       %initialise variables/constants
clear all
close all
c1=0.04;
c2=-c1;
c3=0.08;
c4=-c3;
p1=[];


p2=1.5;
p3=1.75;
y3=[];
z1=3;

s2=0.005;
s3=s2;
sw3=00.00375;
                                                                       %p1 signal generator
                                                                       %dummy denotes the times the gate is lifted over 6 hours.
dummy=[1:10 20:30 40:45 60:63 90:130 190:200 210:212 215:219 300:350];
for i=1:(60*6)
    p1(i)=2.5;
    if ismember(i,dummy)
        p1(i)=2.5+0.2;
    end
        
end

h1=z1-p1;
%generate states
%x(t+1)=Ax + Bu+d;  y(t)=Cx+D
A= [1+c2 0 0; 0 1+c4 c3; 1 0 0];
C=[0 1 0];
x=[];
y=[];
y(1)=1;
x(1:3,1)=[0.25; 0.25; 1];% INITIALISE
for k=2:(6*60)
  B=[c1*h1(k); 0; 0]; %encode input in B
  d=[normrnd(0,s2);normrnd(0,s3); 0];
  x(1:3,k)=A*x(1:3,k-1)+B+d;
  D=normrnd(p3,sw3);
  y(k)=C*x(1:3,k)+ D;
end



%part c
yhat=[];
xhat=[];
yhat(1)=1;
xhat(1:3,1)=[0; 0; 0];

F=A;
G=[c1;0;0];
H=[0 1 0];
V = [s2^2 0 0; 0 s3^2 0; 0 0 0];                                                %covariance matrix
W = sw3^2;                                                                      % additive white noise
P=eye(3)*1000;
innov=0;
K=[];
Kp=F*P*(H')*inv(H*P*(H')+W);
K(1:3,1)=Kp;
yhat(2)=1;

                                                                                % in filter form
xf=[];
Kf=[];
xf(1:3,1)=[1;1;1];
xf(1:3,2)=[1;1;1];
Pf=eye(3);

Kf(1:3,1)=Pf*H'*(inv(H*Pf*H'+W));
Kf(1:3,2)=Pf*H'*(inv(H*Pf*H'+W));
error=y(1)-H*xf(1:3,1);

for k=2:length(y)
xf(1:3,k)=xf(1:3,k-1)+Kf(1:3,k-1)*error;%+G*h1(k); %x(n|n)
error=y(k)-H*xf(1:3,k-1)-p3;%e(n)
Kf(1:3,k)=Pf*H'*(inv(H*Pf*H'+W));%Kf(n)
Pf=Pf-Kf(1:3,k)*H*Pf; %P(n|n)
xf(1:3,k)=F*xf(1:3,k)+G*h1(k);%Next state
Pf=F*Pf*F'+V; %P(N)
end
figure(4)
plot(H*xf)
hold on
plot(xf(1,:))
hold on
plot(x(1,:))
hold on
plot(x(2,:))
hold on
plot(y-p3)
title('kalman filter')
legend('kalman estimate h3','kalman est h2','actual h2','actual h3','measurement')
ylim([0,1])


figure(5)
plot(abs((xf(1,1:360)-x(1,:)).^2))
hold on
plot(abs((xf(2,1:360)-x(2,:)).^2))
title('Error/cost')
legend('error h2','error h3')
ylim([0,1])


**Bold** and _Italic_ and `Code` text


```
![Filter-states](https://github.com/ElliotjFitz/Folio-1/blob/gh-pages/Filter%20perf.png?raw=true)


![Cost/Error over time](https://github.com/ElliotjFitz/Folio-1/blob/gh-pages/cost.png?raw=true)

What we see is that after the initial errors in estimation, the filter will stabilise to a minimal error variance after the 50 hour mark. 



