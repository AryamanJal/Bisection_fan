#function that checks if B_{F, G} = cone(v- u: v in vert(F), u in vert(G)) for every pair of facets F, G of given polytope P.

use application "polytope";
use strict;
use warnings;
require "./bis_cone.pl"; #ensure that this file ends with a 1;
require "./ray_description.pl"; #ensure that this file ends with a 1;

sub check_eq{
    my $p = shift;
    my $n = $p->N_FACETS; 
    my @s = (0..$n-1);
    my $r =0;
    foreach my $i (@s){
        foreach my $j (@s){
            my $f = facet($p, $i);
            my $g = facet($p, $j);
            my $B = B($f, $g);
            my $b = b($f, $g);
            if(equal_polyhedra($B, $b)== false){
                $r = $r + 1;
            }
        }
    }
    if ($r == 0){
        print("The two descriptions are equal!");
    }
    else {
        print("The two descriptions are unequal at \$r places");        
    }
    
    
}