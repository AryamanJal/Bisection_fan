#functions to generate a bisection cone based on the pdf Inequality_description_of_bisection_cones
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
    my $A_F = -$C_F->FACETS;
    my $C_G = new Cone(RAYS=>$V_G);
    my $A_G = -$C_G->FACETS;
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
    my $C_F = new Cone(INPUT_RAYS=>$V_F);
    my $C_G = new Cone(INPUT_RAYS=>$V_G);
    my $A_F = -$C_F->FACETS;
    my $A_G = -$C_G->FACETS;
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
    my $A_F = -$C_F->FACETS;
    my $A_G = -$C_G->FACETS;
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

#making the cone \mathcal{B}_{F, G}
sub B
{
    my $F = shift; #facet 1
    my $G = shift; #facet 2
    my $M = new Matrix(M($F, $G));
    my $B = new Cone(INEQUALITIES=>$M); 
    return $B;
}
1; #ending file like this so we can include it elsewhere using the require keyword