function [u, M, F, basis, d] = ScalarProjection( basis_family, target_fun, degree, domain )
variate = symvar( target_fun );
basis = PolynomialBasisFunction( basis_family, degree, variate, domain );
M = AssembleGramMatrix( basis, domain );
F = AssembleForceVector( basis, target_fun, domain );
d = M \ F;
u = symfun( transpose( d ) * basis, variate );
end

function M = AssembleGramMatrix( basis, domain )
M = int( basis .* transpose( basis ), domain );
end

function F = AssembleForceVector( basis, target_fun, domain )
F = int( basis * target_fun, domain );
end