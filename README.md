# Visualising 2 and 3-dimensional bisection fans
This repo contains a Jupyter notebook and executable polymake scripts for visualising $2$ and $3$-dimensional bisection fans and computing bisection cones. The corresponding paper is "Polyhedral combinatorics of bisectors" and is available here: https://arxiv.org/abs/2308.14372. 

### Installation
Install polymake on your machine as detailed [here](https://polymake.org/doku.php/download/start). If you're using a Windows machine, the best way to install it is to do so on WSL2. Instructions are [here](https://docs.google.com/document/d/1pJm5Shye_7nwL4tEx695frccYMfbXHYSpKEHJ_HHEt0/edit). To run the Jupyter notebook containing polymake commands, add the polymake kernel to your Jupyter set-up as outlined [here](https://polymake.org/doku.php/user_guide/howto/jupyter).

### Ways of executing the polymake code
Clone the repository using `git clone https://github.com/AryamanJal/Bisection_fan.git`.

(1) Executing from the terminal:  `cd` into `Bisection_fan` and open polymake with `polymake`. Choose your favourite $2$ or $3$-dimensional, centrally symmetric polytope $P$ and save it to some variable `$p`. For example `$p=cube(3);`. Then run
```
script("./write_vertices.pl"); 
write_vert($p);
```

This writes the vertices of $P$ to `vertices.txt`. Now run `script("./plot.pl");`. If everything has been installed correctly, a browser window should open up where the collection of all the bisection cones of $P$ are visualised simultaneously. (It might be the case that the visualisation fails to load in your default browser; if so, copy the URL and paste it into a different browser.) The bisection fan of the cube for example will look like the image below.

![bisfancube](cube.png)

To generate the bisection fan of the Fedorov solids (the images at the end of the paper), run step (1) on the following polytopes, replacing `$p` with `$p_i` for $i=1,2,3$ or $4$ as necessary:

1. Truncated octahedron: `$p_1 = archimedean_solid("truncated_octahedron");`
2. Hexagonal prism: `$p_2 = prism(n_gon(6), -1, 1);`
3. Hexarhombic dodecahedron: `$p_3 = load("hexarhombic_dodecahedron.poly");`
4. Rhombic dodecahedron: `$p_4 = catalan_solid("rhombic_dodecahedron");`

(2) Execute from the Jupyter notebook: Instructions in the `Visualising bisection fans.ipynb` are self-explanatory.

### Bisection cones

The file `rand_cent_sym.pl` generates a random $3$-dimensional centrally symmetric polytope on `$n` (even number) vertices that you can plot the bisection cones of with:
```
$n = 8;
script("./rand_cent_sym.pl"); 
$p = rand_cs_p($n);
script("./plot.pl"); 
```

The cone $\mathcal{B}_{F, G}$ was computed from its inequality description, a derivation of which is done in the `Inequality_description_of_bisection_cones.pdf`

The file `check_equality.pl` verifies the proposition in the paper that states that if $F, G$ are facets of the polytope $P$ then $\mathcal{B}_{F, G} = \text{cone}(\{v  - u : v \in \text{vert}(F), u \in \text{vert}(G)\}):$
```
$p = cube(3);
script("./check_equality.pl"); 
check_eq($p);
```

The file `bis_cone.pl` contains all the functions (annotated) necessary to generate the bisection cone $\mathcal{B}_{F, G}$ given facets $F, G$ of $P.$
