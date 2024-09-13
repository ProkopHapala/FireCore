# Conversion of Fireball to Classes, C++ and OpenCL

## How to approach it

* we need to see how to split it to objects
    * we need to see which global variables can be localized
        * we can use dependency graph which variables are used by which subroutines
            * this can be build automatically using either AO tools or static code analyzer
                * we can check it by `use configurations, only x,y,z`