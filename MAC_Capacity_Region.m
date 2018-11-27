clear
clc
m=10;   % m receive antennas
n1=20;  % n1 transmit antennas for user 1
n2=20;  % n2 transmit antennas for user 2
h1=randn(m,n1);  % channel matrix for user 1
h2=randn(m,n2);  % channel matrix for user 2
z=eye(m,m);      % covariance matrix of AWGN
p1=1000;  % power constraint for user 1
p2=1000;  % power constraint for user 2
u1=0:0.01:1;
u2=zeros(1,length(u1));
r1=zeros(1,length(u1));
r2=zeros(1,length(u1));

 for n=1:length(u1)
    u2(n)=1-u1(n);
    if u1(n)<=u2(n)
cvx_begin   % find optimal covariance matrices s1 and s2 by maximizing u1*r1+u2*r2
    variable s1(n1,n1);
    variable s2(n2,n2);
    a1=-u1(n)*.5*det_rootn(h1*s1*h1'+h2*s2*h2'+z);
    b1=-(u2(n)-u1(n))*.5*det_rootn(h2*s2*h2'+z);
    c1=u2(n)*.5*det_rootn(z); 
    minimize (a1+b1+c1)
    subject to
        trace(s1)<=p1; 
        trace(s2)<=p2;
        s1==semidefinite(n1,n1);
        s2==semidefinite(n2,n2);
cvx_end
    if u1(n)==0
        s1=0;
    else
    end
r1(n)=.5*log(det(h1*s1*h1'+h2*s2*h2'+z)/det(z))-.5*log(det(h2*s2*h2'+z)/det(z));
r2(n)=.5*log(det(h2*s2*h2'+z)/det(z));
else
cvx_begin  % find optimal covariance matrices s1 and s2 by maximizing u1*r1+u2*r2
    variable s1(n1,n1);
    variable s2(n2,n2);
    a2=-u2(n)*.5*det_rootn(h1*s1*h1'+h2*s2*h2'+z);
    b2=-(u1(n)-u2(n))*.5*det_rootn(h1*s1*h1'+z);
    c2=u1(n)*.5*det_rootn(z);
    minimize (a2+b2+c2)
    subject to
        trace(s1)<=p1; 
        trace(s2)<=p2;
        s1==semidefinite(n1,n1);
        s2==semidefinite(n2,n2);
cvx_end
    if u2(n)==0
        s2=0;
    else
    end
r1(n)=.5*log(det(h1*s1*h1'+z)/det(z));
r2(n)=.5*log(det(h1*s1*h1'+h2*s2*h2'+z)/det(z))-.5*log(det(h1*s1*h1'+z)/det(z));
    end
 end
plot(r1,r2)  % plot the capacity region
grid on
xlabel('R1')
ylabel('R2')
title('Capacity region of a two-user vector multiple access channel')
