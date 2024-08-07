function [u, M, F, basis, d] = BernoulliEulerBeamEquation( spline_space, youngs_modulus, moment_inertia, displacement_constraint, slope_constraint, distributed_load, moment_load, shear_load )
variate = symvar( distributed_load );
SplineBasis = PiecewisePolynomialBasisFunction( spline_space, variate );
basis =  @(deriv) diff( SplineBasis.Basis, variate, deriv );
domain = spline_space.getSplineDomain();
[M, F] = ApplyGoverningEquation( basis, domain, youngs_modulus, moment_inertia, distributed_load );
[M, F] = ApplyBoundaryConditions( M, F, basis, domain, youngs_modulus, moment_inertia, displacement_constraint, slope_constraint, moment_load, shear_load );
d = M \ F;
u = symfun( transpose( d ) * basis(0), variate );
u = simplify( u, Steps=10 );

    function [M, F] = ApplyGoverningEquation( basis, domain, youngs_modulus, moment_inertia, distributed_load )
        M = int( ( (youngs_modulus * moment_inertia ) * basis(4) ) .* ( (youngs_modulus * moment_inertia ) * transpose( basis(4) ) ), domain );
        F = -int( ( (youngs_modulus * moment_inertia ) * basis(4) ) * distributed_load, domain );
    end

    function [M, F] = ApplyBoundaryConditions( M, F, basis, domain, youngs_modulus, moment_inertia, displacement_constraint, slope_constraint, moment_load, shear_load )
        [MD, FD] = DisplacementBoundaryCondition( basis, domain, displacement_constraint );
        [MN, FN] = SlopeBoundaryCondition( basis, domain, slope_constraint );
        [MM, FM] = MomentBoundaryCondition(basis, domain, youngs_modulus, moment_inertia, moment_load );
        [MS, FS] = ShearBoundaryCondition( basis, domain, youngs_modulus, moment_inertia, shear_load );
        M = M + MD + MN + MM + MS;
        F = F + FD + FN + FM + FS;
    end

    function [M, F] = DisplacementBoundaryCondition( basis, domain, displacement_constraint )
        M = subs( basis(0) .* transpose( basis(0) ),   symvar( basis(0) ), domain(1) );
        F = subs( basis(0) .* displacement_constraint, symvar( basis(0) ), domain(1) );
    end

    function [M, F] = SlopeBoundaryCondition( basis, domain, slope_constraint )
        M = limit( basis(1), variate, domain(1), "right" ) .* transpose( limit( basis(1), variate, domain(1), "right" ) );
        F = limit( basis(1), variate, domain(1), "right" ) .* slope_constraint;
    end

    function [M, F] = MomentBoundaryCondition( basis, domain, youngs_modulus, moment_inertia, moment_load )
        M = ( limit( basis(2), variate, domain(2), "left" ) * (youngs_modulus * moment_inertia ) ) .* transpose( limit( basis(2), variate, domain(2), "left" ) );
        F = limit( basis(2), variate, domain(2), "left" ) .* moment_load;
    end

    function [M, F] = ShearBoundaryCondition( basis, domain, youngs_modulus, moment_inertia, shear_load )
        M = ( limit( basis(3), variate, domain(2), "left" ) * (youngs_modulus * moment_inertia ) ) .* transpose( limit( basis(3), variate, domain(2), "left" ) );
        F = -limit( basis(3), variate, domain(2), "left" ) .* shear_load;
    end
end