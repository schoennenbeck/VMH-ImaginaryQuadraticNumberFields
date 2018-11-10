# VMH
Magma and GAP methods for working with arithmetic groups over imaginary quadratic number fields. These methods are described in the following articles:
* [Perfect lattices over imaginary quadratic number fields](http://de.arxiv.org/abs/1304.0559) (O. Braun and R. Coulangeon, Mathematics of Computation, 2015)
* [Computing in arithmetic groups with Voronoi's algorithm](https://arxiv.org/abs/1407.6234) (O. Braun, R. Coulangeon, G. Nebe and S. Schönnenbeck, Journal of Algebra, 2015)
* [Resolutions for unit groups of orders](https://arxiv.org/abs/1609.08835) (S. Schönnenbeck, Journal of Homotopy and related structures, 2016)

The implementation was primarily done by Oliver Braun and Sebastian Schönnenbeck during their time as doctoral students at RWTH Aachen University funded by the DFG Research Training Group "Experimental and constructive algebra".

Some results of computations performed using these algorithms are already available in our [database](http://www.math.rwth-aachen.de/~Oliver.Braun/unitgroups/).

## Overview
The algorithms we implemented deal with groups of the form $GL_n(K)$ for $n=1,2$ and $K$ an imaginary number field. The following algorithms are available:
* Computing a presentation (i.e. generators and defining relations).
* Solving constructive membership problems (i.e. writing a given element in these generators).
* Computing a contractible G-complex which can be used for homology computations with the GAP-package [HAP](http://hamilton.nuigalway.ie/Hap/www/). These algorithms are also available for finite index subgroups of the base group (as long as we know the index and how to check for membership in this group).

## General Workflow
Overall, the usability of the algorithms is not great and a lot of the general design of the Magma-package is pretty terrible. However, for reasons of me no longer working in academia and generally not wanting to do the refactoring it is going to stay this way for the foresesable future. If you are willing to overlook these weaknesses, performing computations with our packages works as follows:
1. Decide on the field you want to work over and the dimension you have in mind (due to the computational complexity only 2 and 3 actually work).
2. Edit the file BasicData.m accordingly (i.e. provide the dimension n, and the negative integer d whose square-root you want to have in the field)
3. Start Magma and attach the package.
4. Start by calling V:=VoronoiAlgorithm(); This is the initial computation that is needed for everything that follows.
5. If you are only interested in presentations and constructive membership you can do all your computations in Magma. Otherwise you will use some functions (described later) to write some files describing a certain combinatorial structure.
6. These files can be read using the GAP-functions we provide and can be used to obtain (or write to file) contractible G-complexes for use with HAP.
