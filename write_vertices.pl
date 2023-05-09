#Write the vertices of polytope P to the vertices.txt file. 
use application "polytope";
use strict;
use warnings;

$home = $ENV{HOME};

sub write_vert{
    my $p = shift;
    open(my $f, ">", "$home/Polymake_4.9/polymake_projects/Bisection_fan/vertices.txt");
    print $f $p->VERTICES;
    close($f);
}


