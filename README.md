# sphep
`sphep` is a tiny C++ tool for generating evenly distributed points on a d-dimensional sphere. It works by generating random points on the sphere and redistributing them by a constrained n-body simulation.
Requirements:
- Eigen3 library
### test
The generator can be tested using the spheptest inside of the test subfolder.
The output points can be visualized using `gnuplot`:
```splot "points.txt" pointtype 7```
### license
The code (except `FindEigen3.cmake`) is licensed under the MIT-license. See the LICENSE file. The cmake module `FindEigen3.cmake` is licensed under BSD 2-clause license (see the top of that file)
