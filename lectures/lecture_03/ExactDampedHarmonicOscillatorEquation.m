function exactSolution = ExactDampedHarmonicOscillatorEquation( mass, damping, stiffness, force, displacement_constraint, slope_constraint, domain )
syms y(t)
eqn = mass*diff( y, t, 2 ) + damping*diff( y, t, 1 ) + stiffness*diff( y, t, 0 ) == force;
D1y = diff( y, t, 1 );
cond = [     y( domain(1) ) == displacement_constraint, ...
    D1y( domain(1) ) == slope_constraint ];
U(t) = dsolve( eqn, cond );

exactSolution.f = force;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = domain;
exactSolution.U = simplify( U, Steps=10 );
end