use application "polytope";
use strict;
use warnings;
my $home = $ENV{HOME};
#require "/Users/aryamanjal/Polymake_4.6/Bisection_fan/bis_cone.pl"; #ensure that that files end with a 1;
require "$home/Polymake_4.9/polymake_projects/Bisection_fan/bis_cone.pl";

#read the contents of the vertices.txt file and store it to matrix $B.
open(INPUT, "<", "$home/Polymake_4.9/polymake_projects/Bisection_fan/vertices.txt");
my $b = new Matrix<Rational>(<INPUT>);
close(INPUT);

my $p = new Polytope(POINTS=>$b);

#function that returns intersection of \mathcal{B}_{F, G} and polytope P
sub B1
{   
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $M = new Matrix(M($F, $G));
    my $B = new Cone(INEQUALITIES=>(zero_vector($M->rows)|$M)); #we do this so $B and $p have the same ambient dim, a necessary step for the next line.
    #print($B->CONE_AMBIENT_DIM, "\n");
    #print($p->CONE_AMBIENT_DIM);
    my $I = intersection($B, $p);
    return $I,
}

#function that plots the collection of bisection cones.
sub vis
{
    my $n = $p->N_FACETS; 
    my @s = (0..$n-1);
    my @G;
    foreach my $i (@s){
        foreach my $j (@s){
            my $f = facet($p, $i);
            my $g = facet($p, $j);
            my $I = B1($f, $g);
            push(@G, $I);          
            }
        }
    return @G;
}


my @L = vis();
#my @L1 = map($_->VISUAL(VertexLabels=>'hidden'),@L);
#compose(@L1);
#print("hello");
# if($p->DIM == 2){
#         foreach my $i (@s){
#             foreach my $j (@s){
#                 my $f = facet($p, $i);
#                 my $g = facet($p, $j);
#                 my $B = B($f, $g);
#                 print($B->CONE_AMBIENT_DIM, "\n");
#                 print($p->CONE_AMBIENT_DIM);
#                 my $I = intersection($B, $p);
#                 push(@G, $I);          
#             }
#         }
#     return @G;
#     }
compose(map($_->VISUAL(VertexLabels=>'hidden'),@L));

1; #ending file like this so we can include it elsewhere using the require keyword