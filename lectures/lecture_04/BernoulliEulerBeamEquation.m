function [u, M, F, basis, d] = BernoulliEulerBeamEquation( basis, youngs_modulus, moment_inertia, displacement_constraint, slope_constraint, distributed_load, moment_load, shear_load, boundary_tol, penalty )
variate = symvar( distributed_load );
domain = basis.splineSpace.getSplineDomain();
basis = @(deriv) diff( basis.Basis, variate, deriv );
[M, F] = ApplyGoverningEquation( basis, domain, youngs_modulus, moment_inertia, distributed_load );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, youngs_modulus, moment_inertia, displacement_constraint, slope_constraint, moment_load, shear_load, boundary_tol, penalty );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
end

function [M, F] = ApplyGoverningEquation( basis, domain, youngs_modulus, moment_inertia, distributed_load )
M = int( basis(4) * (youngs_modulus * moment_inertia ) .* transpose( basis(4) ), domain );
F = -int( basis(4) * distributed_load, domain );
end

function [M, F] = ApplyBoundaryConditions( M, F, basis, domain, youngs_modulus, moment_inertia, displacement_constraint, slope_constraint, moment_load, shear_load, boundary_tol, penalty )
[MD, FD] = DisplacementBoundaryCondition( basis, domain, displacement_constraint, boundary_tol, penalty );
[MN, FN] = SlopeBoundaryCondition( basis, domain, slope_constraint, boundary_tol, penalty );
[MM, FM] = MomentBoundaryCondition(basis, domain, youngs_modulus, moment_inertia, moment_load, penalty, boundary_tol );
[MS, FS] = ShearBoundaryCondition( basis, domain, youngs_modulus, moment_inertia, shear_load, penalty, boundary_tol );
M = M + MD + MN + MM + MS;
F = F + FD + FN + FM + FS;
end

function [M, F] = DisplacementBoundaryCondition( basis, domain, displacement_constraint, boundary_tol, penalty )
M = penalty * int( basis(0) .* transpose( basis(0) ),   [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = penalty * int( basis(0) .* displacement_constraint, [domain(1) - boundary_tol, domain(1) + boundary_tol] );
end

function [M, F] = SlopeBoundaryCondition( basis, domain, slope_constraint, boundary_tol, penalty )
M = penalty * int( basis(1) .* transpose( basis(1) ), [domain(1) - boundary_tol, domain(1) + boundary_tol] );
F = penalty * int( basis(1) .* slope_constraint,      [domain(1) - boundary_tol, domain(1) + boundary_tol] );
end

function [M, F] = MomentBoundaryCondition( basis, domain, youngs_modulus, moment_inertia, moment_load, penalty, boundary_tol )
M = penalty * int( basis(2) * (youngs_modulus * moment_inertia ) .* transpose( basis(2) ), [domain(2) - boundary_tol, domain(2) + boundary_tol] );
F = penalty * int( basis(2) .* moment_load,                                                [domain(2) - boundary_tol, domain(2) + boundary_tol] );
end

function [M, F] = ShearBoundaryCondition( basis, domain, youngs_modulus, moment_inertia, shear_load, penalty, boundary_tol )
M = penalty * int( basis(3) * (youngs_modulus * moment_inertia ) .* transpose( basis(3) ), [domain(2) - boundary_tol, domain(2) + boundary_tol] );
F = -penalty * int( basis(3) .* shear_load,                                                [domain(2) - boundary_tol, domain(2) + boundary_tol] );
end