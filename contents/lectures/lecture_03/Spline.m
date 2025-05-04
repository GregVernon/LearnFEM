classdef Spline
    % Base properties
    properties
        spline_space SplineSpace
    end

    % Derived properties
    properties
        local_dc_extraction_operators
        local_c0_extraction_operators
        spline_dc_extraction_operator
        spline_c0_extraction_operator
        dc_basis
        basis
    end

    % Constructors
    methods
        function obj = Spline( spline_space )
            obj.spline_space = spline_space;
            obj = obj.InitializeDCBasis();
            obj.basis = obj.ComputeSplineBasis();
        end
    end

    % Initialization
    methods
        function obj = InitializeDCBasis( obj )
            variate = obj.GetVariate();
            basis_name = obj.GetBasisName();
            num_elems = obj.GetNumElements();
            num_dc_basis = obj.ComputeNumDCBasis();
            obj.dc_basis = sym( zeros( num_dc_basis, 1 ) );
            for elem_id = 1 : num_elems
                elem_domain = obj.GetElementDomain( elem_id );
                elem_degree = obj.GetElementDegree( elem_id );
                local_dc_basis_ids = obj.GetSupportedDCBasisIdsFromElementId( elem_id );
                local_basis = PolynomialBasisFunction( basis_name, elem_degree, variate, elem_domain );
                for n = 1 : length( local_basis )
                    if elem_id == num_elems
                        obj.dc_basis(local_dc_basis_ids(n)) = piecewise( elem_domain(1) <= variate & variate <= elem_domain(2), local_basis(n), sym( 0 ) );
                    else
                        obj.dc_basis(local_dc_basis_ids(n)) = piecewise( elem_domain(1) <= variate & variate < elem_domain(2), local_basis(n), sym( 0 ) );
                    end
                end
            end
        end
    end

    % Get information about the spline space
    methods
        function basis_name = GetBasisName( obj )
            basis_name = obj.spline_space.GetBasisName();
        end

        function variate = GetVariate( obj )
            variate = obj.spline_space.GetVariate();
        end

        function vertex = GetVertex( obj, vertex_id )
            vertex = obj.spline_space.GetVertex( vertex_id );
        end

        function num_elems = GetNumElements( obj )
            num_elems = obj.spline_space.GetNumElements();
        end

        function num_interfaces = GetNumInterfaces( obj )
            num_interfaces = obj.spline_space.GetNumInterfaces();
        end

        function domain = GetSplineDomain( obj )
            domain = obj.spline_space.GetSplineDomain();
        end
    end

    % Get information about individual elements of the spline
    methods
        function elem_domain = GetElementDomain( obj, elem_id )
            elem_domain = obj.spline_space.GetElementDomain( elem_id );
        end

        function elem_degree = GetElementDegree( obj, elem_id )
            elem_degree = obj.spline_space.GetElementDegree( elem_id );
        end

        function elem_continuity = GetElementContinuity( obj, elem_id )
            elem_continuity = obj.spline_space.GetElementContinuity( elem_id );
        end

        function elem_interface_ids = GetElementInterfaceIds( obj, elem_id )
            elem_interface_ids = obj.spline_space.GetElementInterfaceIds( elem_id );
        end

        function supported_dc_basis_ids = GetSupportedDCBasisIdsFromElementId( obj, elem_id )
            start_id = sum( obj.spline_space.degree(1:elem_id-1) + 1 ) + 1;
            num_local_dc_basis = obj.GetElementDegree( elem_id ) + 1;
            supported_dc_basis_ids = start_id : start_id + ( num_local_dc_basis - 1 );
        end

        function supported_basis_ids = GetSupportedBasisIdsFromElementId( obj, elem_id )
            supported_basis_ids = [];
            for e = 1 : elem_id
                elem_degree = GetElementDegree( obj, e );
                elem_interface_ids = GetElementInterfaceIds( obj, e );
                left_elem_continuity = GetInterfaceContinuity( obj, elem_interface_ids(1) );
                if e == 1
                    curr_basis_id = 0;
                else
                    curr_basis_id = curr_basis_id - ( left_elem_continuity + 1 );
                end
                for n = 1 : elem_degree + 1
                    curr_basis_id = curr_basis_id + 1;
                    if e == elem_id
                        supported_basis_ids = [ supported_basis_ids; curr_basis_id ];
                    end
                end
            end
        end
    end

    % Get information about individual interfaces of the spline
    methods
        function interface_elems = GetInterfaceElementIds( obj, interface_id )
            interface_elems = obj.spline_space.GetInterfaceElementIds( interface_id );
        end

        function interface_cont = GetInterfaceContinuity( obj, interface_id )
            interface_cont = obj.spline_space.GetInterfaceContinuity( interface_id );
        end

        function interface_vertex_id = GetInterfaceVertexId( obj, interface_id )
            interface_vertex_id = obj.spline_space.GetInterfaceVertexId( interface_id );
        end
    end

    % Compute information about the spline space
    methods
        function knot_vector = ComputeKnotVector( obj )
            knot_vector = obj.spline_space.ComputeKnotVector();
        end

        function num_dc_basis = ComputeNumDCBasis( obj )
            num_dc_basis = obj.spline_space.ComputeNumDCBasis();
        end
    end

    % USpline
    methods

        function spline_basis = ComputeSplineBasis( obj )
            spline_assembly_operator = obj.ComputeSplineAssemblyOperator();
            spline_basis = spline_assembly_operator * obj.dc_basis;
        end

        function spline_assembly_operator = ComputeSplineAssemblyOperator( obj )
            if obj.GetNumElements == 1
                num_dc_basis = obj.ComputeNumDCBasis();
                spline_assembly_operator = sym( eye( num_dc_basis ) );
            else
                global_constraint_matrix = obj.ComputeGlobalInterfaceConstraintMatrix();
                contiguous_index_sets = obj.ComputeGlobalContiguousIndexSets();
                for ii = 1 : length( contiguous_index_sets )
                    A = global_constraint_matrix( :, contiguous_index_sets{ii} );
                    add_vec = obj.generate_independent_unit_vector( A );
                    A = [add_vec; A];
                    b = [1; zeros( size( A, 1) - 1, 1)];
                    c = A \ b;
                    N(ii, contiguous_index_sets{ii} ) = c;
                end
                normalize_cols = ( N * transpose( N ) ) \ ( N * ones( size( global_constraint_matrix, 2 ), 1 ) );
                spline_assembly_operator = normalize_cols .* N;
            end
        end

        function interface_index_sets = ComputeLocalInterfaceIndexSets( obj, interface_id )
            num_interfaces = obj.GetNumInterfaces();
            interface_cont = obj.GetInterfaceContinuity( interface_id );
            interface_elems = obj.GetInterfaceElementIds( interface_id );
            index_set_length = interface_cont + 2;
            if interface_id == 1
                interface_r_elem_id = interface_elems;
                r_elem_dc_basis_ids = obj.GetSupportedDCBasisIdsFromElementId( interface_r_elem_id );
                interface_index_sets = r_elem_dc_basis_ids(1);
            elseif interface_id == num_interfaces
                interface_l_elem_id = interface_elems;
                l_elem_dc_basis_ids = obj.GetSupportedDCBasisIdsFromElementId( interface_l_elem_id );
                interface_index_sets = l_elem_dc_basis_ids(end);
            else
                interface_l_elem_id = interface_elems(1);
                l_elem_dc_basis_ids = obj.GetSupportedDCBasisIdsFromElementId( interface_l_elem_id );
                interface_index_sets = zeros( interface_cont + 1, index_set_length );
                offset = -1*interface_cont;
                for set_id = 1 : interface_cont + 1
                    start_id = l_elem_dc_basis_ids(end) + offset;
                    end_id = start_id + (  index_set_length - 1 );
                    interface_index_sets(set_id,:) = start_id : end_id;
                    offset = offset + 1;
                end
            end
        end

        function contiguous_index_sets = ComputeGlobalContiguousIndexSets( obj )
            num_interfaces = obj.GetNumInterfaces();
            if num_interfaces == 2
                candidate_contiguous_index_sets = obj.ComputeLocalInterfaceIndexSets( 1 );
                contiguous_index_sets = {};
                for ii = 1 : size( candidate_contiguous_index_sets, 1 )
                    contiguous_index_sets{ii,1} = candidate_contiguous_index_sets(ii,:);
                end
                return
            end

            index_sets = {};
            for interface_id = 1 : num_interfaces
                interface_index_sets{interface_id} = obj.ComputeLocalInterfaceIndexSets( interface_id );
                for ii = 1 : size( interface_index_sets{interface_id}, 1 )
                    if isempty( index_sets )
                        index_sets = { interface_index_sets{interface_id}(ii,:) };
                    else
                        index_sets = [ index_sets; interface_index_sets{interface_id}(ii,:) ];
                    end
                end
            end

            num_dc_basis = ComputeNumDCBasis( obj );
            for ii = 1 : num_dc_basis
                supported_basis = false;
                for jj = 1 : length( index_sets )
                    if ismember( ii, index_sets{jj} )
                        supported_basis = true;
                    end
                end
                if supported_basis == false
                    index_sets{end+1} = ii;
                end
            end

            candidate_contiguous_index_sets = {};
            for ii = 1 : length( index_sets )
                candidate_contiguous_index_sets{ii} = index_sets{ii};
                for jj = ii + 1 : length( index_sets )
                    if ( index_sets{jj}(1) ) == ( candidate_contiguous_index_sets{ii}(end) )
                        candidate_contiguous_index_sets{ii} = unique( [candidate_contiguous_index_sets{ii}, index_sets{jj}] );
                    end
                end
            end

            contiguous_index_sets = {};
            for ii = 1 : length( candidate_contiguous_index_sets )
                keep_set = true;
                for jj = 1 : length( candidate_contiguous_index_sets )
                    if ii == jj
                        continue
                    else
                        if all( ismember( candidate_contiguous_index_sets{ii}, candidate_contiguous_index_sets{jj} ) )
                            keep_set = false;
                        end
                    end
                end
                if keep_set == true
                    contiguous_index_sets{end+1,1} = candidate_contiguous_index_sets{ii};
                end
            end

            first_idx = zeros( size( contiguous_index_sets ) );
            for ii = 1 : length( contiguous_index_sets )
                first_idx(ii) = contiguous_index_sets{ii}(1);
            end
            [~,sidx] = sort( first_idx );
            contiguous_index_sets = contiguous_index_sets( sidx );
        end

        function constraint_matrix = ComputeLocalInterfaceConstraintMatrix( obj, interface_id )
            variate = obj.GetVariate();
            interface_vertex_id = obj.GetInterfaceVertexId( interface_id );
            interface_vertex = obj.GetVertex( interface_vertex_id );
            interface_cont = obj.GetInterfaceContinuity( interface_id );
            interface_elems = obj.GetInterfaceElementIds( interface_id );
            num_interfaces = obj.GetNumInterfaces();
            if interface_id == 1 || interface_id == num_interfaces
                constraint_matrix = [];
            else
                interface_l_elem_id = interface_elems(1);
                interface_r_elem_id = interface_elems(2);
                l_elem_dc_basis_ids = obj.GetSupportedDCBasisIdsFromElementId( interface_l_elem_id );
                r_elem_dc_basis_ids = obj.GetSupportedDCBasisIdsFromElementId( interface_r_elem_id );
                num_dc_basis = obj.ComputeNumDCBasis();
                constraint_matrix = sym( zeros( interface_cont + 1, num_dc_basis ) );
                for deriv = 0 : interface_cont
                    for term = 0 : interface_cont
                        l_dc_basis_id = l_elem_dc_basis_ids(end - term);
                        r_dc_basis_id = r_elem_dc_basis_ids(1 + term);
                        l_dc_basis = obj.dc_basis(l_dc_basis_id);
                        r_dc_basis = obj.dc_basis(r_dc_basis_id);
                        l_dc_basis_deriv = diff( l_dc_basis, deriv );
                        r_dc_basis_deriv = diff( r_dc_basis, deriv );
                        l_dc_basis_deriv_const = limit( l_dc_basis_deriv, variate, interface_vertex, "left" );
                        r_dc_basis_deriv_const = limit( r_dc_basis_deriv, variate, interface_vertex, "right" );
                        constraint_matrix(deriv+1, l_dc_basis_id) = l_dc_basis_deriv_const;
                        constraint_matrix(deriv+1, r_dc_basis_id) = -1 * r_dc_basis_deriv_const;
                    end
                end
            end
        end

        function global_constraint_matrix = ComputeGlobalInterfaceConstraintMatrix( obj )
            num_interfaces = obj.GetNumInterfaces();
            for interface_id = 1 : num_interfaces
                local_constraint_matrix = obj.ComputeLocalInterfaceConstraintMatrix( interface_id );
                if interface_id == 1
                    global_constraint_matrix = local_constraint_matrix;
                else
                    global_constraint_matrix = [ global_constraint_matrix; local_constraint_matrix ];
                end
            end
        end

        function vec = generate_independent_unit_vector( obj, matrix )
            [~, num_cols] = size( matrix );
            init_rank_null = rank( null( matrix ) );
            success = false;
            for col = 1 : num_cols
                vec = sym( zeros( 1, num_cols ) );
                vec(col) = 1;
                test_matrix = [ vec; matrix ];
                if rank( null( test_matrix ) ) < init_rank_null
                    success = true;
                    break
                end
            end
            if success == false
                vec(:) = 0;
            end
        end
    end

end