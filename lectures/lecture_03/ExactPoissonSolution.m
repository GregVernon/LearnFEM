function exactSolution = ExactPoissonSolution( f, g, h, domain )
syms u(x)
eqn = diff( u, x, 2 ) + f == 0;
Du = diff( u, x, 1 );
cond = [  u( domain(2) ) == g, ...
         Du( domain(1) ) == h ];
U(x) = dsolve( eqn, cond );

exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = domain;
exactSolution.U = simplify( U, Steps=10 );
end