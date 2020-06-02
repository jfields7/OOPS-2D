# OOPS-2D
The two-dimensional Object-Oriented PDE Solver (OOPS-2D) is an extension of the [original OOPS code](https://github.com/jfields7/OOPS). It's designed to be an easy-to-use scientific code that is written in modern, flexible C++ to help alleviate the terrible tendency of scientists to write hard-coded, difficult-to-read programs that gradually turn into horrible monstrosities from generations of research students modifying and extending the code to solve new problems.

To be fair, we've all done this. We start a research problem, realize that the code can't do what we want it to do, and throw in a few kluges and hacks to make it work well enough for our problem. Unfortunately, this means the code is usually difficult to read, hasn't been integrated well with features that weren't used for your specific problem, and is missing documentation.

By providing a base core of components that all interact seamlessly and creating a flexible interface for your own projects, OOPS-2D makes it a lot easier to work with your own problems. It's not a flawless solution; the code is currently missing a lot of features that are necessary for production codes, such as alternative coordinate systems and grid layouts, mesh refinement, and a complex suite of data analysis tools. The code is 100% geared toward initial boundary value problems using classic Runge-Kutta or multistage solvers, so your mileage may vary with complex adaptive solvers, symplectic integrators, or certain implicit solvers. However, it's the perfect solution for rapid prototyping or learning the basics of scientific computing.

# Documentation
For the OOPS-2D documentation, please see the [wiki](https://github.com/jfields7/OOPS/wiki).
