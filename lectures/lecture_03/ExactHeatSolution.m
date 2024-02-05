function exactSolution = ExactHeatSolution( internal_heat_flux, prescribed_temperature, boundary_heat_flux, domain )
syms u(x)
eqn = diff( u, x, 2 ) + internal_heat_flux == 0;
Du = diff( u, x, 1 );
cond = [  u( domain(2) ) == prescribed_temperature, ... 
         Du( domain(1) ) == boundary_heat_flux ];
U(x) = dsolve( eqn, cond );
 
exactSolution.f = internal_heat_flux;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = domain;
exactSolution.U = U;
end
