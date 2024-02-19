function l2_error = ComputeL2Error( fun1, fun2, domain )
    l2_error = sqrt( int( ( fun1 - fun2 )^2, domain ) );
end