classdef SplineSpace
    % Base properties
    properties
        variate
        vertices
        degree
        continuity
        basis_name
    end
    
    % Constructor
    methods
        function obj = SplineSpace( basis_name, variate, degree, vertices, continuity )
            obj.basis_name = basis_name;
            obj.variate = variate;
            obj.degree = degree;
            obj.vertices = vertices;
            obj.continuity = continuity;
        end        
    end

    % Get information about the spline space
    methods
        function basis_name = GetBasisName( obj )
            basis_name = obj.basis_name;
        end

        function variate = GetVariate( obj )
            variate = obj.variate;
        end
        
        function vertex = GetVertex( obj, vertex_id )
            vertex = obj.vertices( vertex_id );
        end

        function num_elems = GetNumElements( obj )
            num_elems = length( obj.degree );
        end
        
        function num_interfaces = GetNumInterfaces( obj )
            num_interfaces = length( obj.continuity );
        end

        function domain = GetSplineDomain( obj )
            domain = [ obj.vertices(1), obj.vertices(end) ];
        end
    end

    % Get information about individual elements of the spline
    methods
        function elem_domain = GetElementDomain( obj, elem_id )
            elem_domain = obj.vertices( elem_id : elem_id + 1 );
        end
        
        function elem_degree = GetElementDegree( obj, elem_id )
            elem_degree = obj.degree( elem_id );
        end

        function elem_continuity = GetElementContinuity( obj, elem_id )
            elem_continuity = obj.continuity( elem_id : elem_id + 1 );
        end
    end
    
    % Get information about individual interfaces of the spline
    methods
        function interface_elems = GetInterfaceElementIds( obj, interface_id )
            num_interfaces = obj.GetNumInterfaces();
            if interface_id == 1
                interface_elems = interface_id;
            elseif interface_id == num_interfaces
                interface_elems = interface_id - 1;
            else
                interface_elems = [ interface_id - 1, interface_id ];
            end
        end

        function interface_cont = GetInterfaceContinuity( obj, interface_id )
            interface_cont = obj.continuity( interface_id );
        end

        function interface_vertex_id = GetInterfaceVertexId( obj, interface_id )
            interface_vertex_id = interface_id;
        end
    end

    % Compute information about the spline space
    methods
        function knot_vector = ComputeKnotVector( obj )
            knot_vector = cell( 1,length(obj.vertices) );
            for ii = 1:length(obj.vertices)
                knot_vector{ii} = repmat( obj.vertices(ii), 1, obj.degree(1) - obj.continuity(ii) );
            end
            knot_vector = cell2mat( knot_vector );
        end

        function num_dc_basis = ComputeNumDCBasis( obj )
            num_dc_basis = sum( obj.degree + 1 );
        end
    end
end