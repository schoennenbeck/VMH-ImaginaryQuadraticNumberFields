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
