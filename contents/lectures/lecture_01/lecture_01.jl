module lecture_01

    # IMPORTS
    import LinearAlgebra
    import Symbolics
    import SymbolicNumericIntegration
    import Integrals
    import QuadGK
    integrate = QuadGK.quadgk
    import Plots

    include( "lecture_01_utils.jl" )
    using .lecture_01_utils

    # EXPORTS
    export ChangeOfVariable
    export ScalarProjection
    export ComputeL2Error
    export PolynomialBasis
    export PolynomialBasisFunction
    export PlotPolynomialBasis

    # GLOBALS
    @enum PolynomialBasis begin
        Bernstein
        Chebyshev
        Lagrange
        Legendre
        Monomial
    end

    BASIS_ID_STRINGS = Dict( Bernstein => "B",
                             Chebyshev => "T",
                             Lagrange  => "L",
                             Legendre  => "P",
                             Monomial  => "M" )

    # FUNCTION DEFINITIONS

    ## Change of Variable
    function ChangeOfVariable( from_value, from_domain, to_domain )
        ## SHIFT TO ZERO
        to_value = from_value - from_domain[1]
        ## SCALE TO EQUAL RANGE
        to_value *= ( to_domain[2] - to_domain[1] ) / ( from_domain[2] - from_domain[1] )
        ## APPLY TO DOMAIN
        to_value += to_domain[1]
        return to_value
    end

    ## Basis Functions
    function PolynomialBasisFunction( basis_name, degree, variate, domain )
        if basis_name == Bernstein
            basis = BernsteinBasis( degree, variate, domain )
        elseif basis_name == Chebyshev
            basis = ChebyshevBasis( degree, variate, domain )
        elseif basis_name == Lagrange
            basis = LagrangeBasis( degree, variate, domain )
        elseif basis_name == Legendre
            basis = LegendreBasis( degree, variate, domain )
        elseif basis_name == Monomial
            basis = MonomialBasis( degree, variate, domain )
        end
        return basis
    end

    function BernsteinBasis( degree, variate, domain )
        variate = ChangeOfVariable( variate, domain, [0, 1] )
        basis = []
        for i = 0 : degree
            bin_coeff = binomial( degree, i )
            monomial = variate^i
            polynomial = ( 1 - variate ) ^ ( degree - i )
            append!( basis, bin_coeff * monomial * polynomial )
        end
        return basis
    end

    function ChebyshevBasis( degree, variate, domain )
        variate = ChangeOfVariable( variate, domain, [-1, 1] )
        basis = []
        for i = 0 : degree
            if i == 0
                append!( basis, variate^0 )
            elseif i == 1
                append!( basis, variate^1 )
            else
                append!( basis, ( 2 * variate * basis[i] ) - basis[i-1] )
            end
        end
        return basis
    end

    function LagrangeBasis( degree, variate, domain )
        node = LinRange( domain[1], domain[2], degree + 1 )
        basis = []
        for i = 1 : degree + 1
            append!( basis, variate^0 )
            for j = 1 : degree + 1
                if i != j
                    basis[i] *= ( ( variate - node[j] ) / ( node[i] - node[j] ) )
                end
            end
        end
        return basis
    end

    function LegendreBasis( degree, variate, domain )
        variate = ChangeOfVariable( variate, domain, [-1, 1] )
        basis = []
        for i = 0 : degree
            if i == 0
                append!( basis, variate^0 )
            elseif i == 1
                append!( basis, variate^1 )
            else
                n = i - 1
                term_1 = ( 2 * n + 1 ) * variate * basis[i]
                term_2 = n * basis[i-1]
                append!( basis, ( term_1 - term_2 ) / ( n + 1 ) )
            end
        end
        return basis
    end

    function MonomialBasis( degree, variate, domain )
        variate = ChangeOfVariable( variate, domain, [0, 1] )
        basis = []
        for n in range( 0, degree )
            append!( basis, variate^n )
        end
        return basis
    end
    
    ## SCALAR PROJECTION
    function ScalarProjection( target_fun, basis_name, degree, domain )

        function AssembleGramMatrix()
            integrand = basis * transpose( basis )
            integrand = Symbolics.build_function( integrand, variate, expression=false )
            D = integrate.( integrand, domain[1], domain[2]; atol=1e-16, order=12 )
            return D
        end

        function AssembleForceVector()
            integrand = basis * target_fun
            integrand = Symbolics.build_function( integrand, variate, expression=false )
            F = integrate.( integrand, domain[1], domain[2]; atol=1e-16, order=12 )
            return F
        end

        variate = Symbolics.Num( Symbolics.get_variables( target_fun )[1] )
        basis = PolynomialBasisFunction( basis_name, degree, variate, domain )
        basis_dim = length( basis )
        M = AssembleGramMatrix()
        F = AssembleForceVector()
        M = to_numeric( Float64, M )
        F = to_numeric( Float64, F )
        d = vec( M \ F )
        u = transpose( d ) * basis
        return u, M, F, basis, d
    end

    function ComputeL2Error( f1, f2, domain )
        error_fun = (f1 - f2 )^2
        variate = Symbolics.Num( Symbolics.get_variables( error_fun )[1] )
        error_fun_compiled = Symbolics.build_function( error_fun, variate, expression=false )
        l2_error = QuadGK.quadgk( error_fun_compiled, domain[1], domain[2]; atol=1e-16, order=12 )[1]
        return l2_error
    end

    ## PLOTTING
    function PlotPolynomialBasis( basis_name, degree, variate, domain::Vector{Symbolics.Num} )
        basis = PolynomialBasisFunction( basis_name, degree, variate, domain )
        basis_id_string = BASIS_ID_STRINGS[basis_name]
        domain = to_numeric( Float64, domain )
        plot = Plots.plot( basis[1], domain[1], domain[2], label="\$ $(basis_id_string)_0 \$" )
        for i = 1 : degree
            n = i + 1
            Plots.plot!( plot, basis[n], domain[1], domain[2], label="\$ $(basis_id_string)_$(i) \$" )
        end
        return plot
    end

    function PlotPolynomialBasis( basis_name, degree, variate, domain )
        basis = PolynomialBasisFunction( basis_name, degree, variate, domain )
        basis_id_string = BASIS_ID_STRINGS[basis_name]
        plot = Plots.plot( basis[1], domain[1], domain[2], label="\$ $(basis_id_string)_0 \$" )
        for i = 1 : degree
            n = i + 1
            Plots.plot!( plot, basis[n], domain[1], domain[2], label="\$ $(basis_id_string)_$(i) \$" )
        end
        return plot
    end
end