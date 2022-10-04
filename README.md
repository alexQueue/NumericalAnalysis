# Numerical Analysis

By Alex Quinlan, Henry Jacobson, Pia Callmer and Nicola Sabbadini

[![Tests](https://github.com/alexQueue/NumericalAnalysis/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/alexQueue/NumericalAnalysis/actions/workflows/test.yml)

![Latest beam animation](/img/framework_dynamic.gif)

## Running the code
Run the dynamic case from the root of the directory with
```bash
julia --project=. src/single_dynamic.jl
```
```bash
julia --project=. src/framework_dynamic.jl
```
Run the static code with
```bash
julia = --project=. src/single_static.jl
```
```bash
julia = --project=. src/dynamic_static.jl
```
## Testing

Tests are run in Github Actions but are non-blocking. You can run tests locally with
```bash
julia --project=. test/runtests.jl
````

## Package management

See: https://pkgdocs.julialang.org/v1/environments/

Add packages by:
* Entering the Julia shell.
* Then enter package mode with `]`.
* Activate the local project with `activate .`.
* Finally, add the package with `add <package-name>`.
