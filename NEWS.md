v0.10.2
- adding fixed mean model
- adding data set: threeGaussians
- fixing compilation issues on some platforms

v0.10.1
------
- rewriting all C the code in C++11
- adding split method
- adding threads

v0.9.4
------
- adding README.md
- lots of refactoring
- small fixes

v0.9.3
------
- giving up support of -1 iterations (fixing memcheck problems)
- changing the way initial centers vector is handled: for each start, length(centers) clusterings are performed
- adding two datasets: fourGaussians and mixShapes

v0.9.2
------
- checking input data for NA values (session crushing)
- changing 'ZERO_EPSILON' to '1.0e-32'
- lots of refactoring

v0.9.1
------
Initial Release.
