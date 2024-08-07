function [u, A, F, basis, d] = DampedHarmonicOscillatorEquation( basis_family, degree, domain, mass, damping, stiffness, force, displacement_constraint, slope_constraint )
variate = symvar( force );
if basis_family == "Trig"
basis = @(deriv) diff( TrigonometricBasis( degree, variate, domain ), variate, deriv );
else
basis = @(deriv) diff( PolynomialBasisFunction( basis_family, degree, variate, domain ), variate, deriv );
end
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
        A = subs( basis(1) .* transpose( basis(1) ), symvar( basis(0) ), domain(1) );
        F = subs( basis(1) .* slope_constraint,      symvar( basis(0) ), domain(1) );
    end

end