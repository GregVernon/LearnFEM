function [u, M, F, basis, d] = ThirdOrderPDE( basis, distributed_load, displacement_constraint, slope_constraint, gradient_constraint, boundary_tol, penalty )
variate = symvar( distributed_load );
domain = basis.splineSpace.getSplineDomain();
basis = @(deriv) diff( basis.Basis, variate, deriv );
M = AssembleLinearOperatorMatrix( basis, domain );
F = AssembleForceVector( basis, distributed_load, domain );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, displacement_constraint, slope_constraint, gradient_constraint, boundary_tol, penalty );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
end

function M = AssembleLinearOperatorMatrix( basis, domain )
M = int( basis(3) .* transpose( basis(3) ), domain );
end

function F = AssembleForceVector( basis, distributed_load, domain )
F = -1 * int( basis(3) * distributed_load, domain );
end

function [M, F] = ApplyBoundaryConditions( M, F, basis, domain, displacement_constraint, slope_constraint, gradient_constraint, boundary_tol, penalty )
[MD, FD] = DisplacementBoundaryCondition( basis, domain, displacement_constraint, boundary_tol, penalty );
[MN, FN] = SlopeBoundaryCondition( basis, domain, slope_constraint, boundary_tol, penalty );
[MG, FG] = GradientBoundaryCondition(basis, domain, gradient_constraint, boundary_tol, penalty );
M = M + MD + MN + MG;
F = F + FD + FN + FG;
end

function [M, F] = DisplacementBoundaryCondition( basis, domain, displacement_constraint, boundary_tol, penalty )
M = penalty * int( basis(0) .* transpose( basis(0) ),          [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = penalty * int( basis(0) .* displacement_constraint, [domain(1) - boundary_tol, domain(1) + boundary_tol] );
end

function [M, F] = SlopeBoundaryCondition( basis, domain, slope_constraint, boundary_tol, penalty )
M = penalty * int( basis(1) .* transpose( basis(1) ), [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = penalty * int( basis(1) .* slope_constraint, [domain(1) - boundary_tol, domain(1) + boundary_tol] );
end

function [M, F] = GradientBoundaryCondition( basis, domain, gradient_constraint, boundary_tol, penalty )
M = penalty * int( basis(2) .* transpose( basis(2) ), [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = penalty * int( basis(2) .* gradient_constraint, [domain(1) - boundary_tol, domain(1) + boundary_tol] );
end