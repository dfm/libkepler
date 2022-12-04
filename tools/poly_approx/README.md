Tools for generating polynomial approximations to functions using the [Remez
exchange algorithm](https://en.wikipedia.org/wiki/Remez_algorithm). We use this
here to design efficient algorithms for evaluating trig functions and Kepler
starters. I should note that this isn't designed to be very robust or
general-purpose, but it does the trick for what we need here!

Boost has a nice discussion and implementation of this method:

- https://www.boost.org/doc/libs/1_47_0/boost/math/tools/remez.hpp
- https://www.boost.org/doc/libs/1_80_0/libs/math/doc/html/math_toolkit/remez.html
