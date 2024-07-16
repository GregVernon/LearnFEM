function exactSolution = ExactBernoulliEulerBeamSolution( E, I, f, M, Q, domain )
syms u(x)
eqn = E * I * diff( u, x, 4 ) + f == 0;
D1u = diff( u, x, 1 );
D2u = diff( u, x, 2 );
D3u = diff( u, x, 3 );
cond = [   u( domain(1) ) == 0, ...
         D1u( domain(1) ) == 0, ...
         D2u( domain(2) ) == M / ( E * I ), ...
         D3u( domain(2) ) == -Q / ( E * I ) ];
U(x) = dsolve( eqn, cond );

exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = domain;
exactSolution.U = simplify( U, Steps=10 );
end