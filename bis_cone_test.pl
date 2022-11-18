#functions to generate a bisection cone based on the pdf Inequalities_for_bisection_cone
use application "polytope";
use strict;
use warnings;
 
#calculate distance function given facet F of P 
sub dist 
{
    my $F = shift;
    my $V = $F->VERTICES;
    #remove the leading column of ones
    $V = $V->minor(All, ~[0]);
    #take the transpose so that the i^th column of V is the i^th vertex of F
    my $VT = transpose($V);
    my $z_one = ones_vector($F->N_VERTICES);
    $z_one = new Matrix($z_one);
    my $r = common::solve_left($VT, $z_one);
    return $r;
}

#making matrix A
sub A 
{
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $V_F = $F->VERTICES;
    my $V_G = $G->VERTICES;
    $V_F = $V_F->minor(All, ~[0]);
    $V_G = $V_G->minor(All, ~[0]);
    my $C_F = new Cone(RAYS=>$V_F);
    my $A_F = $C_F->FACETS;
    my $C_G = new Cone(RAYS=>$V_G);
    my $A_G = $C_G->FACETS;
    my $z_1 = dist($F);
    my $z_2 = dist($G);
    my $A = new Matrix($A_F/$A_G/($z_1 -$z_2)/($z_2-$z_1));
    return $A;
}       

# #making cone C

sub C
{ 
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $V_F = $F->VERTICES;
    $V_F = $V_F->minor(All, ~[0]); #discard first column of ones
    my $V_G = $G->VERTICES;
    $V_G = $V_G->minor(All, ~[0]); #discard first column of ones
    my $C_F = new Cone(RAYS=>$V_F);
    my $C_G = new Cone(RAYS=>$V_G);
    my $A_F = $C_F->FACETS;
    my $A_G = $C_G->FACETS;
    my $m1 = $A_F->rows;
    my $m2 = $A_G->rows;  
    my $r = $m1+$m2 +2;
    my $i = unit_matrix($r);
    my $i2 = dense($i); 
    my $A = A($F, $G);
    my $AT = transpose($A);    
    my $X = new Matrix($i2/$AT/-$AT);
    my $C = new Cone(INEQUALITIES=>$X);
    return $C;
}    

#extract rays of C
sub C_R
{
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $C = C($F, $G);
    my $C_R = $C->RAYS;
    return $C_R;  
}

#matrix of inequalities for \mathcal{B}_{F, G}

sub M
{
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $C_R = C_R($F, $G);
    my $V_F = $F->VERTICES;
    $V_F = $V_F->minor(All, ~[0]); #discard first column of ones
    my $V_G = $G->VERTICES;
    $V_G = $V_G->minor(All, ~[0]); #discard first column of ones
    my $C_F = new Cone(RAYS=>$V_F);
    my $C_G = new Cone(RAYS=>$V_G);
    my $A_F = $C_F->FACETS;
    my $A_G = $C_G->FACETS;
    my $A_GT = transpose($A_G);
    my @U = [];
    my @Y = [];
    my @W = [];
    my @X = [];
    my @H = [];
    my $x = $A_F->rows;
    my $y = $A_G->rows;
    my $z_1 = dist($F);
    my $z_2 = dist($G);
    for (my $k = 0; $k< $C_R->rows; $k++){
        $U[$k] = new Matrix($C_R->row($k)->slice(range(0, $x -1))); 
        $Y[$k] = new Matrix($C_R->row($k)->slice(range($x, $x+$y -1)));
        $W[$k] = new Matrix($C_R->row($k)->slice(range($x +$y, $x+$y)));
        $X[$k] = new Matrix($C_R->row($k)->slice(range($x+$y+1, $x+$y+1)));
    }
    for (my $k = 0; $k< $C_R->rows; $k++){          
        $H[$k] = ($A_GT)*(transpose($Y[$k])) - $W[$k]->elem(0, 0)*transpose($z_2) + $X[$k]->elem(0, 0)*transpose($z_2);
    }
    my $M = new Matrix(transpose($H[0]));
    
    for (my $k = 1; $k< $C_R->rows; $k++){            
        $M = $M/transpose($H[$k]);           
    }
    $M = remove_zero_rows($M);
    return $M;
}

sub B
{
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $M = new Matrix(M($F, $G));
    my $B = new Cone(INEQUALITIES=>$M);
}
#function to build matrix A
#     my $F = shift; #facet 1
#     my $G = shift; #facet 2
#     
#     $V_F = $V_F->minor(All, ~[0]); #discard first column of ones
#     my $V_G = $G->VERTICES;
#     $V_G = $V_G->minor(All, ~[0]); #discard first column of ones
#     my $C_F = new Cone(RAYS=>$V_F);
#     my $A_F = $C_F->FACETS; #maybe this needs a minus sign?
#     my $C_G = new Cone(RAYS=>$V_G);
#     my $A_G = $C_G->FACETS;
#     my $m1 = $A_F->rows;
#     my $m2 = $A_G->rows;  


# for (my $i =0; $i <$n; $i ++){
#     my $z_one = ones_vector($ff->cone($i)->N_RAYS);
#     $z_one = new Matrix($z_one);
#     @A[$i] = -$ff->cone($i)->FACETS;
#     @C1[$i] = new Cone(INEQUALITIES=>zero_vector(($ff->cone($i)->FACETS)->rows)|$ff->cone($i)->FACETS);
#     my $M1 = new Matrix($M->row($i));
#     my $M2 = new Matrix(-$M->row($i));
#     $M1 = $M1/$M2;
#     @H1[$i] = new Polytope(INEQUALITIES=>$M1);
#     @p1[$i] = intersection($C1[$i], $H1[$i]);
#     @V[$i]= $p1[$i]->VERTICES;
#     $V[$i] = $V[$i]->minor(All, ~[0]);
#     @VT[$i]= transpose($V[$i]);
#     @z[$i] = solve_left($VT[$i], $z_one);
#     }

# my @D = []; 
# #build matrix A, calling it D here.
# for (my $i =0; $i <$n; $i ++){
#     for (my $j =0; $j <$n; $j ++){
# #         $D[$i][$j] = new Matrix<Rational>($A[$i]);
#         $D[$i][$j] = new Matrix($A[$i]/$A[$j]/($z[$i]- $z[$j])/($z[$j]-$z[$i]));
        
#     }
# }

# #taking transpose
# my @DT = [];
# for (my $i =0; $i <$n; $i ++){
#     for (my $j =0; $j <$n; $j ++){
#         $DT[$i][$j] = transpose($D[$i][$j]);
        
#     }
# }

#matrix of inequalities P for cone C

# my @P = [];
# for (my $i =0; $i <$n; $i ++){
#     for (my $j =0; $j <$n; $j ++){
#         my $x = $A[$i]->rows;
#         my $y = $A[$j]->rows;
#         my $r = $x+$y +2;
#         my $m1 = unit_matrix($r);
#         my $m2 = dense($m1);    
#         $P[$i][$j] = new Matrix($m2/$DT[$i][$j]/-$DT[$i][$j]);        
#     }
# }

# # making cone C

# my @C = [];
# for (my $i =0; $i <$n; $i ++){
#     for (my $j =0; $j <$n; $j ++){
#         $C[$i][$j] = new Cone(INEQUALITIES=>$P[$i][$j]);
        
#     }
# }

# #extracting rays of C
# my @CR = [];
# for (my $i =0; $i <$n; $i ++){
#     for (my $j =0; $j <$n; $j ++){
#         $CR[$i][$j] = new Matrix($C[$i][$j]->RAYS);
        
#     }
# }

# #making matrix of inequalities M(i, j) for cone A_{F_{i}, F_{j}}
# my @U = [];
# my @Y = [];
# my @W = [];
# my @X = [];
# my @H = [];

# for (my $i =0; $i <$n; $i ++){    
#     for (my $j =0; $j <$n; $j ++){
#         my $x = $A[$i]->rows;
#         my $y = $A[$j]->rows;
#         for (my $k = 0; $k< $CR[$i][$j]->rows; $k++){
#             $U[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range(0, $x -1))); 
#             $Y[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range($x, $x+$y -1)));
#             $W[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range($x +$y, $x+$y)));
#             $X[$i][$j][$k] = new Matrix($CR[$i][$j]->row($k)->slice(range($x+$y+1, $x+$y+1)));
#         }
#     }
# }
# for (my $i =0; $i <$n; $i ++){    
#     for (my $j =0; $j <$n; $j ++){
#         for (my $k = 0; $k< $CR[$i][$j]->rows; $k++){        
#             $H[$i][$j][$k] = $Y[$i][$j][$k]*$A[$j] - $W[$i][$j][$k]->elem(0, 0)*$z[$j] + $X[$i][$j][$k]->elem(0, 0)*$z[$j];
#         }
#     }
    
# }

# my @M = [];
# for (my $i =0; $i <$n; $i ++){    
#     for (my $j =0; $j <$n; $j ++){
#         $M[$i][$j] = new Matrix($H[$i][$j][0]); 
#         for (my $k = 1; $k< $CR[$i][$j]->rows; $k++){            
#             $M[$i][$j] = $M[$i][$j]/$H[$i][$j][$k];           
#         }
#         $M[$i][$j] = remove_zero_rows($M[$i][$j]);
#     }
# }

# # print($M[1][2]);

# my @I = [];
# my @s = (0..$n-1);
# foreach my $i (@s){
#     foreach my $j (@s){ 
# #         next if $i==$j;
#         $I[$i][$j] = new Cone(INEQUALITIES=>zero_vector($M[$i][$j]->rows)|$M[$i][$j]);             
        
#     }
# }

# # my @G= ();
# # my @s1 = (0..$n-1);
# # foreach my $i (@s1){
# #     my @s = (0..$n-1);
# #     splice(@s, $i, 1);
# #     foreach my $j (@s1){
# #         next if $i==$j;
# #         push(@G, $I[$i][$j])
# #     }
    
# # }

# my @G= ();
# my @s1 = (0..$n-1);
# foreach my $i (@s1){    
#     foreach my $j (@s1){
# #         next if $i==$j;
#         push(@G, intersection($I[$i][$j], $p))
#     }
    
# }


# # my @p1 = [];
# # for (my $i =0; $i <$F3->rows; $i ++){
# #     my $M1 = new Matrix($F3->row($i));
# #     my $M2 = new Matrix(-$F3->row($i));
# #     $M1 = $M1/$M2;
# #     $p1[$i] = new Polytope(INEQUALITIES=>$M1);  
# #     }


# # push(@G, @p1);
# # compose(map{$_->VISUAL(VertexLabels=>'hidden')}@p1);

# #print join "  ",@G;
# # my @F= ($I[0][1]);
# # @arr= (2..7)
# # foreach my $i (@arr){
# #     push(@F, $I[0,$i]);
# # }

# # compose($I[0][0]->VISUAL, $I[1][1]->VISUAL, $I[2][2]->VISUAL, $I[3][3]->VISUAL, $I[4][4]->VISUAL, $I[5][5]->VISUAL, $I[6][6]->VISUAL, $I[7][7]->VISUAL);

# # print($I[1][2]->AMBIENT_DIM, "\n", cross(3)->AMBIENT_DIM);
# # $pp= intersection(cross(3), $I[1][2]);

# # compose(map{$_->VISUAL(VertexLabels=>'hidden')}@G);

# intersection($I[1][2], $p)->VISUAL;
# print($I[1][2]->RAYS);

# # tikz(compose(map{$_->VISUAL(VertexLabels=>'hidden')}@G), File=>"bisfan"); uncomment for tikz code

# # compose($I[0][1]->VISUAL, $I[7][0]->VISUAL);
# # print scalar(@G);



# # sleep 2;
# # $checkedfan = fan::check_fan_objects(@G);
# # print $checkedfan->MAXIMAL_CONES;

# # print $checkedfan->MAXIMAL_CONES;
# # for (my $i =0; $i < $n; $i ++){
# #     print($z[$i]);

# # }
