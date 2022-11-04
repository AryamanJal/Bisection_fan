# Visualising 2 and 3-dimensional bisection fans
This repo contains a Jupyter notebook and execute polymake code (the .pl files) for visualising 2 and 3-dimensional bisection fans. 

### Installation
Install polymake on your machine as detailed [here](https://polymake.org/doku.php/download/start). If you're using a Windows machine, the best way to install it is to do so on WSL2. Instructions are [here](https://docs.google.com/document/d/1pJm5Shye_7nwL4tEx695frccYMfbXHYSpKEHJ_HHEt0/edit). To be able to run Jupyter notebook containing polymake commands, add the polymake kernel to your Jupyter set-up as outlined [here](https://polymake.org/doku.php/user_guide/howto/jupyter).

### Ways of executing the polymake code
First clone the repository into a directory called `bis_fan`. You then have two choices:
1) Executing from the terminal:  `cd` into `bis_fan` and open polymake with `polymake`. Choose your favourite $2$ or $3$-dimensional, centrally symmetric polytope $P$ and save it to some variable `$p`. Then type 
```
open(my $f, ">", "vertices.txt"); print $f $p->VERTICES; close($f);

```

This writes the vertices of $P$ to `vertices.txt`. Now run `script("bis_cone2.pl")`. If everything has been installed correctly, a browser window should open up containing the collection of bisection cones all visualised at once. The bisection fan of the octahedron for example will look like the image below.

2) Execute from the Jupyter notebook: Instructions in the `Visualising bisection fan.ipynb` are self-explanatory.

