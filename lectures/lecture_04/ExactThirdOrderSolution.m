function exactSolution = ExactThirdOrderSolution( f, a, b, c, domain )
syms u(x)
eqn = diff( u, x, 3 ) + f == 0;
D1u = diff( u, x, 1 );
D2u = diff( u, x, 2 );
cond = [   u( domain(1) ) == a, ... 
         D1u( domain(1) ) == b, ...
         D2u( domain(1) ) == c ];
U(x) = dsolve( eqn, cond );
 
exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = domain;
exactSolution.U = simplify( U, Steps=10 );

end