# wagglefit 0.2.0

This version adds optimisation via c++.

## New features
---

* Scout and scout recruit superposition models are now numerically optimised in Rcpp (C/C++) using the `nlopt` package.

## Other features
---

* Testing has been extended to all functions, however, there is currently an issue with testing the c functions using nlopt. Tests work locally but not on check, failing with `nlopt_create not provided by nloptr` suggesting the linkers are not being built on the check. Despite the test failing the methods work as expected. Current nlopt tests are coded up but not run. Objective is to get this working.

# wagglefit 0.1.0

This is the first versioning of the waggle fit package and contains some functions to model honeybee waggle dance patterns, however this is very limited and acts as an initial sandbox whilst developing.
