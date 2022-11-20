#function that checks if B_{F, G} = cone(v- u: v in vert(F), u in vert(G)) for every pair of facets F, G.

use application "polytope";
use strict;
use warnings;

sub check_eq{
    my $p = shift;
    my $n = $p->N_FACETS; 
}