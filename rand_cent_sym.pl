#function that takes a positive integer n and outputs a random 3d centrally symmetric polytope on n vertices
use strict;
use warnings;
use application "polytope";

sub rand_cs_p
{
    my $n = shift;
    my $p = rand_sphere(3, $n);
    my $V = $p->VERTICES;
    my $M = rand_vert($V, $n/2);
    my $M1 = $M->minor(All, ~[0]);
    $M1 = ones_vector($M1->rows)|-$M1;
    $M = $M/$M1;
    my $q = new Polytope(VERTICES=>$M);
    return $q;
}

