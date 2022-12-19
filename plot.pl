use application "polytope";
use strict;
use warnings;
my $home = $ENV{HOME};
#require "/Users/aryamanjal/Polymake_4.6/Bisection_fan/bis_cone_test.pl"; #ensure that this files end with a 1;
require "$home/Polymake_4.6/Bisection_fan/bis_cone_test.pl";

#read the contents of the vertices.txt file and store it to matrix $B.
open(INPUT, "<", "$home/Polymake_4.6/Bisection_fan/vertices.txt");
my $B = new Matrix<Rational>(<INPUT>);
close(INPUT);

my $p = new Polytope(POINTS=>$B);

#function that returns intersection of \mathcal{B}_{F, G} and polytope P
sub B1
{
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $M = new Matrix(M($F, $G));
    my $B = new Cone(INEQUALITIES=>zero_vector($M->rows)|$M); 
    my $I = intersection($B, $p);
    return $I,
}

#function that plots the collection of bisection cones.
sub vis
{
    my $n = $p->N_FACETS; 
    my @s = (0..$n-1);
    my @G = ();
    foreach my $i (@s){
        foreach my $j (@s){
            my $f = facet($p, $i);
            my $g = facet($p, $j);
            my $I = B1($f, $g);
            push(@G, $I);          
        }
    }
    compose(map{$_->VISUAL(VertexLabels=>'hidden')}@G);
}

