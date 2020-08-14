## Welcome to GitHub Pages

### Markdown
The following is the implementation of a Kalman Filter coded without using ss and kalmf. 
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
V = [s2^2 0 0; 0 s3^2 0; 0 0 0]; %covariance matrix
W = sw3^2; % additive white noise
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

![Filter_in_action](/Filter%20perf.png)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/ElliotjFitz/github-slideshow/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.

You can use the [editor on GitHub](https://github.com/ElliotjFitz/github-slideshow/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

