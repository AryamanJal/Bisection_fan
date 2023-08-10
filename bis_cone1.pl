use application "polytope";
use strict;
use warnings;

open(INPUT, "<", "/home/aryaman/polymake_notebooks/vertices.txt");
my $B = new Matrix<Rational>(<INPUT>);
close(INPUT);

my $p = new Polytope(POINTS=>$B);

my $F1 = $p->FACETS;
my $F2 = $F1->minor(All, ~[0]);
my $z_vec=zero_vector($F2->rows);
my $F3 = new Matrix($z_vec|$F2);

my $ff = fan::face_fan($p);


my $n = $ff->MAXIMAL_CONES->rows;

my @A = [];
my @V = [];

my @z = [];

for (my $i =0; $i <$n; $i ++){
    my $z_one = ones_vector($ff->cone($i)->N_RAYS);
    $z_one = new Matrix($z_one);
    @A[$i] = -$ff->cone($i)->FACETS;
    @V[$i] = $ff->cone($i)->RAYS;
    @z[$i] = solve_left(transpose($V[$i]), $z_one); 
#     $z_one*inv($V[$i]);
    $z[$i] = new Matrix($z[$i]);

}


my @D = []; 
#build matrix A, calling it D here.
for (my $i =0; $i <$n; $i ++){
    for (my $j =0; $j <$n; $j ++){
#         $D[$i][$j] = new Matrix<Rational>($A[$i]);
        $D[$i][$j] = new Matrix($A[$i]/$A[$j]/($z[$i]- $z[$j])/($z[$j]-$z[$i]));
        
    }
}

#taking transpose
my @DT = [];
for (my $i =0; $i <$n; $i ++){
    for (my $j =0; $j <$n; $j ++){
        $DT[$i][$j] = transpose($D[$i][$j]);
        
    }
}

#matrix of inequalities P for cone C

my @P = [];
for (my $i =0; $i <$n; $i ++){
    my $x = $A[$i]->rows;
    my $r = 2*$x +2;
    my $m1 = unit_matrix($r);
    my $m2 = dense($m1);
    for (my $j =0; $j <$n; $j ++){
        $P[$i][$j] = new Matrix($m2/$DT[$i][$j]/-$DT[$i][$j]);
        
    }
}

# making cone C

my @C = [];
for (my $i =0; $i <$n; $i ++){
    for (my $j =0; $j <$n; $j ++){
        $C[$i][$j] = new Cone(INEQUALITIES=>$P[$i][$j]);
        
    }
}

#extracting rays of C
my @CR = [];
for (my $i =0; $i <$n; $i ++){
    for (my $j =0; $j <$n; $j ++){
        $CR[$i][$j] = new Matrix($C[$i][$j]->RAYS);
        
    }
}

#making matrix of inequalities M(i, j) for cone A_{F_{i}, F_{j}}
my @U = [];
my @Y = [];
my @W = [];
my @X = [];
my @H = [];

for (my $i =0; $i <$n; $i ++){
    my $x = $A[$i]->rows;
    for (my $j =0; $j <$n; $j ++){
        for (my $k = 0; $k< $CR[$i][$j]->rows; $k++){
            $U[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range(0, $x -1))); 
            $Y[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range($x, 2*$x -1)));
            $W[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range(2*$x, 2*$x)));
            $X[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range(2*$x+1, 2*$x+1)));
        }
    }
}
for (my $i =0; $i <$n; $i ++){    
    for (my $j =0; $j <$n; $j ++){
        for (my $k = 0; $k< $CR[$i][$j]->rows; $k++){        
            $H[$i][$j][$k] = $Y[$i][$j][$k]*$A[$j] - $W[$i][$j][$k]->elem(0, 0)*$z[$j] + $X[$i][$j][$k]->elem(0, 0)*$z[$j];
        }
    }
    
}

my @M = [];
for (my $i =0; $i <$n; $i ++){    
    for (my $j =0; $j <$n; $j ++){
        $M[$i][$j] = new Matrix($H[$i][$j][0]); 
        for (my $k = 1; $k< $CR[$i][$j]->rows; $k++){            
            $M[$i][$j] = $M[$i][$j]/$H[$i][$j][$k];           
        }
        $M[$i][$j] = remove_zero_rows($M[$i][$j]);
    }
}

# print($M[1][2]);

my @I = [];
my @s = (0..$n-1);
foreach my $i (@s){
    foreach my $j (@s){ 
        next if $i==$j;
        $I[$i][$j] = new Cone(INEQUALITIES=>zero_vector($M[$i][$j]->rows)|$M[$i][$j]);             
        
    }
}


# my @G= ();
# my @s1 = (0..$n-1);
# foreach my $i (@s1){
#     my @s = (0..$n-1);
#     splice(@s, $i, 1);
#     foreach my $j (@s1){
#         next if $i==$j;
#         push(@G, $I[$i][$j])
#     }
    
# }

my @G= ();
my @s1 = (0..$n-1);
foreach my $i (@s1){    
    foreach my $j (@s1){
        next if $i==$j;
        push(@G, intersection($I[$i][$j], $p))
    }
    
}

# my @p1 = [];
# for (my $i =0; $i <$F3->rows; $i ++){
#     my $M1 = new Matrix($F3->row($i));
#     my $M2 = new Matrix(-$F3->row($i));
#     $M1 = $M1/$M2;
#     $p1[$i] = new Polytope(INEQUALITIES=>$M1);  
#     }


# push(@G, @p1);
# compose(map{$_->VISUAL(VertexLabels=>'hidden')}@p1);

#print join "  ",@G;
# my @F= ($I[0][1]);
# @arr= (2..7)
# foreach my $i (@arr){
#     push(@F, $I[0,$i]);
# }

# compose($I[0][1]->VISUAL, $I[0][2]->VISUAL, $I[0][3]->VISUAL, $I[0][4]->VISUAL, $I[0][5]->VISUAL, $I[0][6]->VISUAL, $I[0][7]->VISUAL);

# print($I[1][2]->AMBIENT_DIM, "\n", cross(3)->AMBIENT_DIM);
# $pp= intersection(cross(3), $I[1][2]);
compose(map{$_->VISUAL(VertexLabels=>'hidden')}@G);

# tikz(compose(map{$_->VISUAL(VertexLabels=>'hidden')}@G), File=>"bisfan");

# compose($I[0][1]->VISUAL, $I[7][0]->VISUAL);
# print scalar(@G);

# sleep 2;
# $checkedfan = fan::check_fan_objects(@G);
# print $checkedfan->MAXIMAL_CONES;

# print $checkedfan->MAXIMAL_CONES;
# for (my $i =0; $i < $n; $i ++){
#     print($z[$i]);

# }
