#!/usr/bin/env julia

using Pkg

Pkg.activate(dirname(@__FILE__))

Pkg.test(; julia_args = ["-O0"]) # Compile time dominates the tests, this helps