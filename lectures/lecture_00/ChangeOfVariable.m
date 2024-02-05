function x = ChangeOfVariable( x, from_domain, to_domain )
x = x - from_domain(1);
x = x * ( ( to_domain(2) - to_domain(1) ) / ( from_domain(2) - from_domain(1) ) );
x = x + to_domain(1);
end