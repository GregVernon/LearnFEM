classdef PiecewisePolynomialBasisFunction

    properties
        splineSpace SplineSpace
        Variate
        Basis
    end

    methods
        function obj = PiecewisePolynomialBasisFunction( splineSpace, variate )
            obj.splineSpace = splineSpace;
            obj.Variate = variate;
            spline_type = obj.splineSpace.getSplineType();
            switch spline_type
                case "BSpline"
                    obj.Basis = obj.buildBsplineBasis();
                case "PiecewiseContinuous"
                    obj.Basis = obj.buildContinuousBasis();
                case "PiecewiseDiscontinuous"
                    obj.Basis = obj.buildDiscontinuousBasis();
            end
        end

        function basis = buildBsplineBasis( obj )
            variate = obj.Variate;
            degree = obj.splineSpace.degree(1);
            kV = obj.splineSpace.getKnotVector();
            N = cell( degree + 1, 1 );
            for p = 0:degree
                N{p+1} = sym( zeros( length( kV ) - ( p + 1 ), 1 ) );
                for ii = 1:length(kV)-(p+1)
                    if p == 0
                        N{p+1}(ii) = piecewise(kV(ii) <= variate & variate < kV(ii+1), sym(1), sym(0) );
                    else
                        divisor(1) = (kV(ii+p)-kV(ii));
                        if divisor(1) == 0
                            term(1) = sym(0);
                        else
                            term(1) = ((variate - kV(ii))/(kV(ii+p)-kV(ii)));
                        end

                        divisor(2) = (kV(ii+p+1)-kV(ii+1));
                        if divisor(2) == 0
                            term(2) = sym(0);
                        else
                            term(2) = ((kV(ii+p+1) - variate)/(kV(ii+p+1)-kV(ii+1)));
                        end
                        N{p+1}(ii) = term(1)*N{p}(ii) + term(2)*N{p}(ii+1);
                    end
                end
            end
            basis = simplify(N{end});
            domain = obj.splineSpace.getSplineDomain();
            basis(end) = piecewise( variate == domain(end), 1, basis(end) );
        end

        function basis = buildContinuousBasis( obj )
            num_elem = obj.splineSpace.getNumElements();
            for e = 1 : num_elem
                elem_degree = obj.splineSpace.degree(e);
                elem_domain = obj.splineSpace.getElementDomain( e );
                elem_basis = PolynomialBasisFunction( "Lagrange", elem_degree, obj.Variate, elem_domain );
                if e == 1
                    curr_basis_ids = ( 1 : elem_degree + 1 );
                    basis = sym( zeros( elem_degree + 1, 1 ) );
                else
                    curr_basis_ids = ( curr_basis_ids(end) : curr_basis_ids(end) + elem_degree );
                    basis( curr_basis_ids(2:end) ) = sym( zeros( elem_degree, 1 ) );
                end
                if e == num_elem
                    condition = ( obj.Variate >= elem_domain(1) ) & ( obj.Variate <= elem_domain(2) );
                else
                    condition = ( obj.Variate >= elem_domain(1) ) & ( obj.Variate < elem_domain(2) );
                end
                for n = 1 : elem_degree + 1
                    basis(curr_basis_ids(n)) = basis(curr_basis_ids(n)) + piecewise( condition, elem_basis(n), sym(0) );
                end
            end
        end

        function basis = buildDiscontinuousBasis( obj )
            num_elem = obj.splineSpace.getNumElements();
            basis = sym( zeros( sum( obj.splineSpace.degree + 1 ), 1 ) );
            for e = 1 : num_elem
                elem_degree = obj.splineSpace.degree(e);
                elem_domain = obj.splineSpace.getElementDomain( e );
                elem_basis = PolynomialBasisFunction( "Chebyshev", elem_degree, obj.Variate, elem_domain );
                if e == 1
                    curr_basis_ids = ( 1 : elem_degree + 1 );
                else
                    curr_basis_ids = ( curr_basis_ids(end) + 1 : curr_basis_ids(end) + elem_degree + 1 );
                end
                if e == num_elem
                    condition = ( obj.Variate >= elem_domain(1) ) & ( obj.Variate <= elem_domain(2) );
                else
                    condition = ( obj.Variate >= elem_domain(1) ) & ( obj.Variate < elem_domain(2) );
                end
                for n = 1 : elem_degree + 1
                    basis(curr_basis_ids(n)) = basis(curr_basis_ids(n)) + piecewise( condition, elem_basis(n), sym(0) );
                end
            end
        end
    end
    
    methods
        function [to_coeffs, R] = SplineChangeOfBasis( obj, to_spline, from_coeffs )
            domain = obj.splineSpace.getSplineDomain();
            D = int( to_spline.Basis .* transpose( to_spline.Basis ), domain );
            C = int( to_spline.Basis .* transpose( obj.Basis ), domain );
            R = D \ C;
            to_coeffs = R * from_coeffs;
        end

        function [to_coeffs, R] = PiecewiseSplinePolynomialChangeOfBasis( obj, basis_name, from_coeffs  )
            num_elem = obj.splineSpace.getNumElements();
            to_coeffs = cell( num_elem, 1 );
            R = cell( num_elem, 1 );
            for e = 1 : num_elem
                degree = obj.splineSpace.degree(e);
                elem_domain = obj.splineSpace.getElementDomain( e );
                spline_basis_ids = obj.basisSupportedByElement( e );
                poly_basis = PolynomialBasisFunction( basis_name, degree, obj.Variate, elem_domain );
                D = int( poly_basis .* transpose( poly_basis ), elem_domain );
                C = int( poly_basis .* transpose( obj.Basis(spline_basis_ids) ), elem_domain );
                R{e} = D \ C;
                to_coeffs{e} = R{e} * from_coeffs(spline_basis_ids);
            end
        end

        function supported_elems = elementSupportedByBasis( obj, basis_id )
            supported_elems = [];
            parts = children( obj.Basis( basis_id ) );
            domains = transpose( [ parts{:, 2} ] );
            for ii = 1 : length( domains )
                domain = children( domains( ii ) );
                if length( domain ) > 1
                    domain = children( domain{ 2 } );
                    domain = [ domain{:} ];
                    if length( domain ) > 1
                        xc = double( sum( domain ) / 2 );
                        supported_elems = [ supported_elems find( obj.splineSpace.vertices <= xc, 1, 'last' ) ];
                    end
                end
            end
        end
        
        function supported_basis = basisSupportedByElement( obj, elem_id )
            supported_basis = [];
            num_basis = length( obj.Basis );
            for ii = 1 : num_basis
                supported_elems = elementSupportedByBasis( obj, ii );
                if ismember( elem_id, supported_elems )
                    supported_basis = [ supported_basis ii ];
                end
            end            
        end
    end
end