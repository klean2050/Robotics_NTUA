function X = trajectory(t0,tf,xA,xB,g)
    
    syms a0 a1 a2 a3 a4 a5 b0 b1 c0 c1 c2 c3 c4 c5
    t1 = t0+(tf-t0)*0.1; t2 = t0+(tf-t0)*0.9;
    eqn1 = a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5 - xA == 0;
    eqn2 = a1 + 2*a2*t0 + 3*a3*t0^2 + 4*a4*t0^3 + 5*a5*t0^4 == 0;
    eqn3 = 2*a2 + 6*a3*t0 + 12*a4*t0^2 + 20*a5*t0^3 - g == 0;
    eqn4 = a0 + a1*t1 + a2*t1^2 + a3*t1^3 + a4*t1^4 + a5*t1^5 - (b0 + b1*t1) == 0;
    eqn5 = a1 + 2*a2*t1 + 3*a3*t1^2 + 4*a4*t1^3 + 5*a5*t1^4 -b1 == 0;
    eqn6 = 2*a2 + 6*a3*t1 + 12*a4*t1^2 + 20*a5*t1^3 == 0;
    eqn7 = c0 + c1*t2 + c2*t2^2 + c3*t2^3 + c4*t2^4 + c5*t2^5 - (b0 + b1*t2) == 0;
    eqn8 = c1 + 2*c2*t2 + 3*c3*t2^2 + 4*c4*t2^3 + 5*c5*t2^4 - b1 == 0;
    eqn9 = 2*c2 + 6*c3*t2 + 12*c4*t2^2 + 20*c5*t2^3 == 0;
    eqn10 = c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4 + c5*tf^5 - xB == 0;
    eqn11 = c1 + 2*c2*tf + 3*c3*tf^2 + 4*c4*tf^3 + 5*c5*tf^4 == 0;
    eqn12 = 2*c2 + 6*c3*tf + 12*c4*tf^2 + 20*c5*tf^3 +g == 0;
    eqn13 = b0+b1*t1 - (a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5) - (0.5*(t1-t0)*b1) == 0;
    eqn14 = (c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4 + c5*tf^5) - (b0+b1*t2) - (0.5*(t1-t0)*b1) == 0;
    
    [A,B] = equationsToMatrix([ eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9 eqn10 eqn11 eqn12 eqn13 eqn14 ] , [ a0 a1 a2 a3 a4 a5 b0 b1 c0 c1 c2 c3 c4 c5 ]);
    X = vpa(linsolve(A,B));
    
end
