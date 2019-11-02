function Jinv = inverse(l0,l1,l2,l3,l4,l5,q1,q2,q3)

    % FORWARD KINEMATIC ANALYSIS
    A00s = [1,0,0,l1 ; 0,0,-1,0 ; 0,1,0,-l0 ; 0,0,0,1];
    A0s1 = [sin(q1),0,cos(q1),l3*sin(q1) ; -cos(q1),0,sin(q1),-l3*cos(q1) ; 0,-1,0,0 ; 0,0,0,1];
    A12 = [cos(q2),-sin(q2),0,l4*cos(q2) ; sin(q2),cos(q2),0,l4*sin(q2) ; 0,0,1,l2 ; 0,0,0,1];
    A2E = [cos(q3),-sin(q3),0,l5*cos(q3) ; sin(q3),cos(q3),0,l5*sin(q3) ; 0,0,1,0 ; 0,0,0,1];
    A01 = A00s*A0s1; A02 = A01*A12; T = A02*A2E;
    % DIFFERENTIAL ANALYSIS
    z0s = A00s(1:3,3); z1 = A01(1:3,3); z2 = A02(1:3,3);
    p0s = A00s(1:3,4); p1 = A01(1:3,4); p2 = A02(1:3,4); pe = T(1:3,4);
    J = [cross(z0s,pe-p0s),cross(z1,pe-p1),cross(z2,pe-p2) ; z0s,z1,z2]; Jinv = inv(J(1:3,:));

end
