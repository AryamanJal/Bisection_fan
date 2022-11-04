# Visualising 2 and 3-dimensional bisection fans
This repo contains a Jupyter notebook and execute polymake code (the .pl files) for visualising 2 and 3-dimensional bisection fans. 

### Installation
Install polymake on your machine as detailed [here](https://polymake.org/doku.php/download/start). If you're using a Windows machine, the best way to install it is to do so on WSL2. Instructions are [here](https://docs.google.com/document/d/1pJm5Shye_7nwL4tEx695frccYMfbXHYSpKEHJ_HHEt0/edit). To be able to run Jupyter notebook containing polymake commands, add the polymake kernel to your Jupyter set-up as outlined [here](https://polymake.org/doku.php/user_guide/howto/jupyter).

### Ways of executing the polymake code
First clone the repository into a directory of your choice. You then have two choices:
1) Executing from the terminal:  `cd` into the directory you cloned into and and open polymake with `polymake`. Run `script(bis_cone2)`
The Visualising bisection fan Jupyter notebook contains some snippets, examples and rough work of polymake code.  The file bis_cone.pl can be run from the command line (after opening polymake) as script("some folder/bis_cone.pl"); where you need to replace some folder with the location where bis_cone.pl was downloaded. e.g. /Dell/Downloads/bis_cone.pl and so on. The bis_cone1.txt contains the same content as bis_cone but is a text file.
