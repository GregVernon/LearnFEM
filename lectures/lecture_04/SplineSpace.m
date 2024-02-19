classdef SplineSpace
    %SPLINESPACE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        degree
        vertices
        continuity
    end

    methods
        function obj = SplineSpace( degree, vertices, continuity )
            %SPLINESPACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.degree = degree;
            obj.vertices = vertices;
            obj.continuity = continuity;
        end
    end

    % Get information about the spline space
    methods
        function num_elems = getNumElements( obj )
            num_elems = length( obj.degree );
        end

        function domain = getSplineDomain( obj )
            domain = [ obj.vertices(1), obj.vertices(end) ];
        end

        function knotVector = getKnotVector( obj )
            knotVector = cell( 1,length(obj.vertices) );
            for ii = 1:length(obj.vertices)
                knotVector{ii} = repmat( obj.vertices(ii), 1, obj.degree(1) - obj.continuity(ii) );
            end
            knotVector = cell2mat( knotVector );
        end

        function spline_type = getSplineType( obj )
            if numel( unique( obj.degree ) ) > 1 && any( obj.continuity > 0 )
                spline_type = "USpline";
            elseif numel( unique( obj.degree ) ) == 1 && any( obj.continuity > 0 )
                spline_type = "BSpline";
            elseif all( obj.continuity == -1 )
                spline_type = "PiecewiseDiscontinuous";
            elseif all( obj.continuity <= 0 )
                spline_type = "PiecewiseContinuous";
            end
        end

        function supported_basis = getSupportedPolynomialBasis( obj )
            spline_type = obj.getSplineType();
            switch spline_type    
                case "USpline"
                    supported_basis = [ "Lagrage", "Bernstein", ]
                case "BSpline"
                case "PiecewiseContinuous"
                case "PiecewiseDiscontinuous"
            end
        end
    end

    % Get information about individual elements of the spline
    methods
        function elem_domain = getElementDomain( obj, elem_id )
            elem_domain = obj.vertices( elem_id : elem_id + 1 );
        end
    end
end