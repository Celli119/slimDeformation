# slimDeformation
This repo is a matlab wrapper for the [slim](http://cs.nyu.edu/~panozzo/papers/SLIM-2016.pdf) algorithm.
The mex function enables to deform a 3d mesh (tetrahedron) without any flipped elements.

![SLIM](http://igl.ethz.ch/publications/2016-slim.jpg "SLIM")


Build instructions:
- `git clone  --recursive https://github.com/Celli119/slimDeformation.git`
- `mkdir build`
- `cd build`
- `cmake ..`
- `make`
