#Write the vertices of polytope P to the vertices.txt file. 
use application "polytope";
use strict;
use warnings;

sub write_vert{
    my $p = shift;
    open(my $f, ">", "./vertices.txt");
    print $f $p->VERTICES;
    close($f);
}


