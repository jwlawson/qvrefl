# qvrefl

Check whether any semi-positive quasi-Cartan companions serve vertices of quivers.

### Background

A quiver *Q* can be represented as a skew-symmetric matrix. A quasi-Cartan
companion of such a quiver matrix *B(Q) = (b<sub>ij</sub>)* is a symmetric matrix *A(Q) =
(a<sub>ij</sub>)* such that diagonal entries are all 2 and off diagonal entries
of *A* have the same absolute values as the matching values in *B,* i.e.
*|a<sub>ij</sub>| = |b<sub>ij</sub>|*.

Any quiver has a number of quasi-Cartan companions, differing by the signs of
the off diagonal values.

A quasi-Cartan matrix *A* gives a quadratic form on some vector space, and the
quiver-mutations can be generalised in some sense to mutations on vectors in
this space. Mutation in the *k*-th direction acts on *(v1, ... , vn)* to give
```
v'i = vi - (vi, vk)vk     if i -> k in Q
      vi                  otherwise.
```
This gives another set of vectors, which have a gram matrix *G = ( (v'i, v'j) )*
using the quadratic form defined by *A*.

The quasi-Cartan matrix *A* **serves** the vertex *k* if *G* is a quasi-Cartan
companion of the mutated matrix *&mu;<sub>k</sub>B = B'*, i.e. 
*|g<sub>ij</sub>| = |b'<sub>ij</sub>|*.

This program finds which, if any of the semi-positive definite (i.e. all
eigenvalues are non-negative) quasi-Cartan matrices of a given quiver serve its
vertices.

### Usage

```
qvrefl -m matrix [-h]
  -m Specify the matrix used to find quasi-Cartan companions
	-h Show this usage information
```

### Output

The program will compute which semi-positive quasi-Cartan companion serves each
vertex of the provided quiver matrix. The output is then the first matrix it
encounters for each vertex, or a statement that no such semi-positive matrix was
found.

Example:
```
>>>qvrefl -m "{ { 0 1 0 0 } { -1 0 1 1 } { 0 -1 0 1 } { 0 -1 -1 0 } }"
0 served by
        2        1        0        0
        1        2        1        1
        0        1        2        1
        0        1        1        2
1 served by
        2        1        0        0
        1        2        1        1
        0        1        2        1
        0        1        1        2
2 is not served by any semi-positive quasi-Cartan
3 served by
        2        1        0        0
        1        2        1        1
        0        1        2        1
        0        1        1        2
```

### Build

Run `make` to build `qvrefl`. To build and run the unit tests use `make test`.

##### Dependencies

The main program depends on the following libraries:

 - [libqv]
 - [armadillo]
 - BLAS/LAPACK

In addition, the tests use the [googletest] framework.

### License

`qvrefl` is available under the Apache v2.0 license.

```
Copyright 2016 John Lawson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```

[libqv]: https://github.com/jwlawson/qv
[armadillo]: http://arma.sourceforge.net/
[googletest]: https://github.com/google/googletest
