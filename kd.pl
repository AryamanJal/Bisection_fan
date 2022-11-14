#functions that takes an integer and returns the symmetric edge polytope corresponding to K_d
use application "polytope";
sub k_d {
    my $n = shift;
    my @r = [];
    my $z_vec = zero_vector($n);
    my $M = new Matrix($z_vec);
    
    for (my $i =0; $i <$n; $i ++){
        for (my $j =0; $j <$n; $j ++){
#             print(unit_vector($n, $i)- unit_vector($n, $j));
            $r[$i][$j] = new Matrix(unit_vector($n, $i)- unit_vector($n, $j));
            $M = new Matrix($M/$r[$i][$j]);            
        }
    }
    $M = remove_zero_rows($M);
    my $z_one = ones_vector($M->rows);
    $z_one = new Matrix($z_one);
    $z_one = transpose($z_one);
    $M = $z_one|$M;
    my $p = new Polytope(VERTICES=>$M);
    return $p;    
}