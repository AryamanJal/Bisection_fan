#function that generates the cone cone(v-u: v in vert(F), u in vert(G))

use application "polytope";
use strict;
use warnings;

sub b{
    my $f = shift;
    my $g = shift;
    my $V_f = $f->VERTICES;
    $V_f = $V_f->minor(All, ~[0]);
    my $V_g = $g->VERTICES;
    $V_g = $V_g->minor(All, ~[0]);
    my $n = $V_f->rows;
    my $m = $V_g->rows;
    my @s1 = (0..$n-1);
    my @s2 = (0..$m-1);
    my @v = [];
    my @u = [];
    my @w = [];

    foreach my $i (@s1){
        foreach my $j (@s2){
            $v[$i] = $V_f->row($i);
            $u[$j] = $V_g->row($j);
            $w[$i][$j] = $v[$i] - $u[$j]; 
            $w[$i][$j] = new Matrix($w[$i][$j]);        
        }
    }
    my $r = new Matrix($w[0][0]);
    foreach my $i (@s1){
        foreach my $j (@s2){
            $r = $r/$w[$i][$j];                    
        }
    }
    $r = $r->minor(~[0], All);
    $r = remove_zero_rows($r);
    # $r = zero_vector($r->rows)|$r;
    # print($w[0][0], "\n");
    # print($r->row(0));

    my $c = new Cone(INPUT_RAYS=>$r);
    return $c;
}

1; #ending file like this so we can include it elsewhere using the require keyword