function [u, M, F, basis, d] = ScalarProjection( basis_family, target_fun, degree, domain )
variate = symvar( target_fun );
basis = PolynomialBasisFunction( basis_family, degree, variate, domain );
M = AssembleGramMatrix( basis, domain );
F = AssembleForceVector( basis, target_fun, domain );
d = M \ F;
u = symfun( transpose( d ) * basis, variate );
end

function M = AssembleGramMatrix( basis, domain )
integrand = simplify( basis .* transpose( basis ), Steps=100 );
M = int( integrand, domain );
M = simplify( M, Steps=100 );
end

function F = AssembleForceVector( basis, target_fun, domain )
integrand = simplify( basis * target_fun, Steps=100 );
F = int( integrand, domain );
F = simplify( F, Steps=100 );
end