clear all; clc; close all; pi = sym(pi);
% *** PART A *** %

%%% FORWARD KINEMATIC ANALYSIS %%%

l0 = sym("l0");
l1 = sym("l1");
l2 = sym("l2");
l3 = 0;
l4 = sym("l4");
l5 = sym("l5");
q1 = sym("q1");
q2 = sym("q2");
q3 = sym("q3");

A00s = [1,0,0,l1 ; 0,0,-1,0 ; 0,1,0,-l0 ; 0,0,0,1];
A0s1 = [sin(q1),0,cos(q1),l3*sin(q1) ; -cos(q1),0,sin(q1),-l3*cos(q1) ; 0,-1,0,0 ; 0,0,0,1];
A12 = [cos(q2),-sin(q2),0,l4*cos(q2) ; sin(q2),cos(q2),0,l4*sin(q2) ; 0,0,1,l2 ; 0,0,0,1];
A2E = [cos(q3),-sin(q3),0,l5*cos(q3) ; sin(q3),cos(q3),0,l5*sin(q3) ; 0,0,1,0 ; 0,0,0,1];

A01 = A00s*A0s1;
A01 = simplify(A01);
A02 = A01*A12;
A02 = simplify(A02);
T = A02*A2E;
T = simplify(T);

%%% DIFFERENTIAL ANALYSIS %%%

z0s = A00s(1:3,3);
z1 = A01(1:3,3);
z2 = A02(1:3,3);
p0s = A00s(1:3,4);
p1 = A01(1:3,4);
p2 = A02(1:3,4);
pe = T(1:3,4);

syms s1 s2 s3 w1 w2 w3 u1 u2 u3 w
% Forward
J = [cross(z0s,pe-p0s),cross(z1,pe-p1),cross(z2,pe-p2) ; z0s,z1,z2];
J = simplify(J);
Jdet = simplify(det(J(1:3,:)));
% Singularities
vars = [q1 q2 q3];
[s1,s2,s3,parameters,conditions] = solve(Jdet==0,vars,"ReturnConditions",true);
Jdet = l4*l5*sin(q3)*(l4*cos(q2)+l5*cos(q2+q3)); % final simplified Jdet
% Inverse
Jinv = simplify(inv(J(1:3,:))); % could also = adjoint(J(1:3,:))/Jdet
w = Jinv*[u1 u2 u3]';           % angular velocities from known end-point velocities
