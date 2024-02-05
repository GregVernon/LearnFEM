function [u, M, F, basis, d] = ScalarProjection( basis_family, target_fun, degree, domain, boundary_tol, penalty )
variate = symvar( target_fun );
basis = PolynomialBasisFunction( basis_family, degree, variate, domain );
M = AssembleGramMatrix( basis, domain, boundary_tol, penalty );
F = AssembleForceVector( basis, target_fun, domain, boundary_tol, penalty );
d = M \ F;
u = symfun( transpose( d ) * basis, variate );
end

function M = AssembleGramMatrix( basis, domain, boundary_tol, penalty )
M = int( basis .* transpose( basis ), [domain(1) + boundary_tol, domain(2) - boundary_tol] );
M = M + penalty * int( basis .* transpose( basis ), [domain(1) - boundary_tol, domain(1) + boundary_tol] );
M = M + penalty * int( basis .* transpose( basis ), [domain(2) - boundary_tol, domain(2) + boundary_tol] );
end

function F = AssembleForceVector( basis, target_fun, domain, boundary_tol, penalty )
F = int( basis * target_fun, [domain(1) + boundary_tol, domain(2) - boundary_tol] );
min_boundary_val = subs( target_fun, symvar( target_fun ), domain(1) );
max_boundary_val = subs( target_fun, symvar( target_fun ), domain(2) );
F = F + penalty * int( basis * min_boundary_val, [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = F + penalty * int( basis * max_boundary_val, [domain(2) - boundary_tol, domain(2) + boundary_tol] );
end