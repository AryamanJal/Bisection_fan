use application "polytope";
use strict;
use warnings;
my $home = $ENV{HOME};
#require "/Users/aryamanjal/Polymake_4.6/Bisection_fan/bis_cone_test.pl"; #ensure that that files end with a 1;
require "$home/Polymake_4.6/Bisection_fan/bis_cone.pl";
require "$home/Polymake_4.6/Bisection_fan/plot.pl";

#read the contents of the vertices.txt file and store it to matrix $B.
open(INPUT, "<", "$home/Polymake_4.6/Bisection_fan/vertices.txt");
my $b = new Matrix<Rational>(<INPUT>);
close(INPUT);

my $p = new Polytope(POINTS=>$b);

#returns collection of hyperplanes passing through the origin and parallel to a facet

sub opp_bis_cones{
    my $n = $p->N_FACETS; 
    my @s = (0..$n-1);
    my @G;
    foreach my $i (@s){
        my $f = facet($p, $i);
        my $g = scale($f, -1);
        my $B = B1($f, $g);
        push(@G, $B);
    }
    return @G;
}

sub diag_bis_cones{
    my $n = $p->N_FACETS; 
    my @s = (0..$n-1);
    my @G;
    foreach my $i (@s){
        my $f = facet($p, $i);
        my $B = B1($f, $f);
        push(@G, $B);
    }
    return @G;
}