clear all; clc; close all; pi = sym(pi);
% *** PART B *** %

%%% Motion Parameters %%%

% length of robot links (cm)
l0=0.0;  l1=0.0;  l3=0.0;
l2=2.0;  l4=2.5;  l5=1.0;
% sampling period & motion period (sec)
sample = 0.01;
duration = 10.0;
t = sample:sample:2*duration;
% initial & final end-point position & acceleration (cm,sec)
xA = l1+l2;	yA = 0;
xB = 1;     yB = -1;
gx = 0.05;   gy = 0.05;

%%% Desired Task-Space Trajectory %%%

xd  = double(trajectory(0,duration,xA,xB,gx));
yd  = double(trajectory(0,duration,yA,yB,gy));
rxd = double(trajectory(duration,2*duration,xB,xA,-gx));
ryd = double(trajectory(duration,2*duration,yB,yA,-gy));

kmax = duration/sample + 1;
x = zeros(length(t),1); y=x; ux=x; uy=x; uz=x; z=x;
for k = 1:2*kmax-1
   time  = (k-1)*sample;
   if (time<duration*0.1)
        x(k)  = polyval(xd(6:-1:1),time);
        ux(k) = polyval(polyder(xd(6:-1:1)),time);
        y(k)  = polyval(yd(6:-1:1),time);
        uy(k) = polyval(polyder(yd(6:-1:1)),time);
   elseif (time<duration*0.9)
        x(k)  = polyval(xd(8:-1:7),time);
        ux(k) = polyval(polyder(xd(8:-1:7)),time);
        y(k)  = polyval(yd(8:-1:7),time);
        uy(k) = polyval(polyder(yd(8:-1:7)),time);
   elseif (time<duration)
        x(k)  = polyval(xd(14:-1:9),time);
        ux(k) = polyval(polyder(xd(14:-1:9)),time);
        y(k)  = polyval(yd(14:-1:9),time);
        uy(k) = polyval(polyder(yd(14:-1:9)),time);
   elseif (time<duration*1.1)
        x(k)  = polyval(rxd(6:-1:1),time);
        ux(k) = polyval(polyder(rxd(6:-1:1)),time);
        y(k)  = polyval(ryd(6:-1:1),time);
        uy(k) = polyval(polyder(ryd(6:-1:1)),time);
   elseif (time<duration*1.9)
        x(k)  = polyval(rxd(8:-1:7),time);
        ux(k) = polyval(polyder(rxd(8:-1:7)),time);
        y(k)  = polyval(ryd(8:-1:7),time);
        uy(k) = polyval(polyder(ryd(8:-1:7)),time);
   else
        x(k)  = polyval(rxd(14:-1:9),time);
        ux(k) = polyval(polyder(rxd(14:-1:9)),time);
        y(k)  = polyval(ryd(14:-1:9),time);
        uy(k) = polyval(polyder(ryd(14:-1:9)),time);
   end
   z(k) = -l3-l4;
end
x  = x(1:2000);   y = y(1:2000);  z = z(1:2000);
ux = ux(1:2000); uy = uy(1:2000);

fig1 = figure;

subplot(2,2,1); 
plot(t,x); 
ylabel('x (cm)'); 
xlabel('time (sec)');

subplot(2,2,2); 
plot(t,y); 
ylabel('y (cm)'); 
xlabel('time (sec)');

subplot(2,2,3); 
plot(t,ux); 
ylabel('ux (cm/sec)'); 
xlabel('time (sec)');

subplot(2,2,4); 
plot(t,uy); 
ylabel('uy (cm/sec)'); 
xlabel('time (sec)');

%%% INVESRE KINEMATICS: DESIRED MOTION -> JOINT SPACE

% for one solution
px2 = sqrt(x.^2+z.^2-l2^2);
q3 = -acos((y.^2+(x-l1).^2+(z+l0).^2-l2^2-l4^2-l5^2)./(2*l4*l5));
q2 = atan2(y,px2)-atan2(l5*sin(q3),(l4+l5*cos(q3)));
q1 = -acos((x*l2-px2.*z)./(px2.^2+l2^2));
qd1 = zeros(length(t),1); qd2=qd1; qd3=qd1;
for k = 1:length(t)
    j = inverse(l0,l1,l2,l3,l4,l5,q1(k),q2(k),q3(k));
    qd1(k) = j(1,1)*ux(k) + j(1,2)*uy(k);
    qd2(k) = j(2,1)*ux(k) + j(2,2)*uy(k);
    qd3(k) = j(3,1)*ux(k) + j(3,2)*uy(k);
end

fig2 = figure;

subplot(2,3,1); 
plot(t,q1);
ylabel('q1 (rad)'); 
xlabel('time (sec)');

subplot(2,3,2); 
plot(t,q2); 
ylabel('q2 (rad)'); 
xlabel('time (sec)');

subplot(2,3,3); 
plot(t,q3);
ylabel('q3 (rad)'); 
xlabel('time (sec)');

subplot(2,3,4); 
plot(t,qd1);
ylabel('qd1 (rad)'); 
xlabel('time (sec)');

subplot(2,3,5); 
plot(t,qd2); 
ylabel('qd2 (rad)'); 
xlabel('time (sec)');

subplot(2,3,6); 
plot(t,qd3);
ylabel('qd3 (rad)'); 
xlabel('time (sec)');

%%% FORWARD KINEMATICS: JOINT MOTION -> CARTESIAN POSITIONS/VELOCITIES

x1 = l2*cos(q1);
y1 = zeros(length(t),1);
z1 = l2.*sin(q1);
x2 = l1 + l2*cos(q1) + l4*cos(q2).*sin(q1);
y2 = l4*sin(q2);
z2 = l2.*sin(q1) - l0 - l4.*cos(q1).*cos(q2);
xe = l1 + l2*cos(q1) + l4*cos(q2).*sin(q1) + l5*cos(q2).*cos(q3).*sin(q1) - l5*sin(q1).*sin(q2).*sin(q3);
ye = l5*sin(q2 + q3) + l4*sin(q2);
ze = l2.*sin(q1) - l0 - l4.*cos(q1).*cos(q2) - l5.*cos(q1).*cos(q2).*cos(q3) + l5.*cos(q1).*sin(q2).*sin(q3);

%%% ANIMATION %%%
 
fig3 = figure;
xlabel('x (cm)');
ylabel('y (cm)');
hold on
plot3([0],[0],[0],'o'); T=5;
for tkn=1:100:T*2*kmax
   tk = mod(tkn,2*kmax-2);
   pause(0.2);
   %hold off    % uncomment this in order to see the actual robot motion
   plot3(0,0,0);
   hold on
   axis([-1 3 -1.5 2.5 -4 1])
   axis on
   plot3(x,y,z,'rs'); 
   plot3([0,x1(tk)],[0,y1(tk)],[0,z1(tk)]);					
   plot3([x1(tk)],[y1(tk)],[z1(tk)],'o');    
   plot3([x1(tk),x2(tk)],[y1(tk),y2(tk)],[z1(tk),z2(tk)]);
   plot3([x2(tk)],[y2(tk)],[z2(tk)],'o'); 
   plot3([x2(tk),xe(tk)],[y2(tk),ye(tk)],[z2(tk),ze(tk)]);
   plot3([xe(tk)],[ye(tk)],[ze(tk)],'y*');
   plot3([x(tk)],[y(tk)],[ze(tk)],'g+');
end
