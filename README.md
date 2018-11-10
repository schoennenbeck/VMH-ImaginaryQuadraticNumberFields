# VMH
Magma and GAP methods for working with arithmetic groups over imaginary quadratic number fields. These methods were primarily developed by Oliver Braun and Sebastian Schönnenbeck during their time as doctoral students at RWTH Aachen University funded by the DFG Research Training Group "Experimental and constructive algebra".

## Overview
The algorithms we implemented deal with groups of the form $GL_n(K)$ for $n=1,2$ and $K$ an imaginary number field. The following algorithms are available:
* Computing a presentation (i.e. generators and defining relations).
* Solving constructive membership problems (i.e. writing a given element in these generators).
* Computing a contractible G-complex which can be used for homology computations with the GAP-package [HAP](http://hamilton.nuigalway.ie/Hap/www/). These algorithms are also available for finite index subgroups of the base group (as long as we know the index and how to check for membership in this group).
