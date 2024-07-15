function [u, M, F, basis, d] = ScalarProjection( basis_family, target_fun, degree, domain )
variate = symvar( target_fun );
basis = PolynomialBasisFunction( basis_family, degree, variate, domain );
M = AssembleGramMatrix( basis, domain );
F = AssembleForceVector( basis, target_fun, domain );
[M, F] = ApplyBoundaryConditions( M, F, basis, target_fun, domain );
d = M \ F;
u = symfun( transpose( d ) * basis, variate );
end

function M = AssembleGramMatrix( basis, domain )
M = int( basis .* transpose( basis ), domain );
end

function F = AssembleForceVector( basis, target_fun, domain )
F = int( basis * target_fun, domain );
end

function [M, F] = ApplyBoundaryConditions( M, F, basis, target_fun, domain )
M = M + subs( basis .* transpose( basis ), symvar( basis ), domain(1) );
M = M + subs( basis .* transpose( basis ), symvar( basis ), domain(2) );
F = F + subs( basis * target_fun, symvar( basis ), domain(1) );
F = F + subs( basis * target_fun, symvar( basis ), domain(2) );
end

