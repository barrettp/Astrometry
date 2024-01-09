# SOFA Routines

The `SOFA` sub-package is meant to be a migration path for those who are already
using the C interface package, i.e., the `sisl` package. This means that `SOFA`
is essentially a drop-in replacement. It also has the potential for checking the
accuracy of the C algorithms by replacing the `Float64` type with higher
precision types such as `Float128` or `BigFloat` and for checking the error
propagation of the algorithms.

In addition to that, it has the potential to provide higher performance compared
to the C library, particularly when doing calculations of a large number (of
order millions) of objects. Having to do a C call on every object incurs a huge
overhead and is a signifficant performance hit.

```@docs
SOFA.Astrom
SOFA.Leapsecond
SOFA.Driftsecond
SOFA.Ldbody
```
