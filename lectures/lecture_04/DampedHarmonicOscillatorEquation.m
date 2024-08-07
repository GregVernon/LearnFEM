function [u, A, F, basis, d] = DampedHarmonicOscillatorEquation( spline_space, mass, damping, stiffness, force, displacement_constraint, slope_constraint )
variate = symvar( force );
SplineBasis = PiecewisePolynomialBasisFunction( spline_space, variate );
basis =  @(deriv) diff( SplineBasis.Basis, variate, deriv );
domain = spline_space.getSplineDomain();
[A, F] = ApplyGoverningEquations( basis, domain, mass, damping, stiffness, force );
[A, F] = ApplyBoundaryConditions( A, F, basis, domain, displacement_constraint, slope_constraint );
d = A \ F;
u = symfun( transpose( d ) * basis(0), variate );
u = simplify( u, Steps=10 );

    function [A, F] = ApplyGoverningEquations( basis, domain, mass, damping, stiffness, force )
        M = int( ( mass * basis(2) ) .* ( mass * transpose( basis(2) ) ), domain );
        C = int( ( mass * basis(2) ) .* ( damping * transpose( basis(1) ) ), domain );
        K = int( ( mass * basis(2) ) .* ( stiffness * transpose( basis(0) ) ), domain );
        A = M + C + K;
        F = int( ( mass * basis(2) ) * force, domain );
    end

    function [A, F] = ApplyBoundaryConditions( A, F, basis, domain, displacement_constraint, slope_constraint )
        [AD, FD] = DisplacementBoundaryCondition( basis, domain, displacement_constraint );
        [AN, FN] = SlopeBoundaryCondition( basis, domain, slope_constraint );
        A = A + AD + AN;
        F = F + FD + FN;
    end

    function [A, F] = DisplacementBoundaryCondition( basis, domain, displacement_constraint )
        A = subs( basis(0) .* transpose( basis(0) ),   symvar( basis(0) ), domain(1) );
        F = subs( basis(0) .* displacement_constraint, symvar( basis(0) ), domain(1) );
    end

    function [A, F] = SlopeBoundaryCondition( basis, domain, slope_constraint )
        A = limit( basis(1), variate, domain(1), "right" ) .* transpose( limit( basis(1), variate, domain(1), "right" ) );
        F = limit( basis(1), variate, domain(1), "right" ) .* slope_constraint;
    end

end