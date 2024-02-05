function x = ChangeOfVariable( x, domain, target_domain )
x = x - domain(1);
x = x * ( ( target_domain(2) - target_domain(1) ) / ( domain(2) - domain(1) ) );
x = x + target_domain(1);
end