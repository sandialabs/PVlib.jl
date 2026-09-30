# PVlib.jl

[![CI](https://github.com/sandialabs/PVlib.jl/actions/workflows/CI.yaml/badge.svg)](https://github.com/sandialabs/PVlib.jl/actions/workflows/CI.yaml)
[![codecov](https://codecov.io/gh/sandialabs/PVlib.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sandialabs/PVlib.jl)

`PVlib.jl` is a Julia package for photovoltaic system modeling and optimization. It provides a Julia-native implementation of selected photovoltaic performance models, with a focus on workflows that benefit from automatic differentiation, gradient-based optimization, and integration with other renewable energy models.

The package includes tools for solar position calculation, plane-of-array irradiance, module temperature, effective irradiance, DC power modeling, and AC power modeling. It also adds utilities for floating photovoltaic applications, including time-varying panel orientation from platform motion, rolling-average cell temperature, ocean surface albedo, and geometric shading from nearby obstacles.

Documentation is available at:

[https://sandialabs.github.io/PVlib.jl/](https://sandialabs.github.io/PVlib.jl/)

## Relationship to pvlib-python

`PVlib.jl` is based in part on [`pvlib-python`](https://pvlib-python.readthedocs.io/), a widely used open-source Python package for photovoltaic system modeling. Several core modeling routines from `pvlib-python` have been adapted and implemented in Julia, including functionality for solar position, irradiance calculations, temperature modeling, and PV power prediction.

`PVlib.jl` is not intended to be a complete Julia replacement for `pvlib-python`. Instead, it provides a focused subset of photovoltaic modeling capabilities designed for Julia-based simulation, automatic differentiation, and gradient-based optimization workflows. This makes it especially useful for applications where PV performance models need to be embedded within larger optimization problems or coupled renewable energy system models.

## Main differences from pvlib-python

The main differences between `PVlib.jl` and `pvlib-python` are:

- **Programming language**:  
  `pvlib-python` is written in Python, while `PVlib.jl` is written in Julia.

- **Scope**:  
  `pvlib-python` is a mature, broad-featured photovoltaic modeling library. `PVlib.jl` currently implements a smaller subset of photovoltaic modeling capabilities, focused on models needed for optimization and floating PV studies. For example, `pvlib-python` includes a broader collection of clear-sky models, irradiance decomposition and transposition methods, single-diode electrical models, tracking-system utilities, loss models, high-level `ModelChain` workflows, and data interfaces.

- **Automatic differentiation**:  
  `PVlib.jl` is designed to work with Julia automatic differentiation tools. Several discontinuous operations used in photovoltaic modeling are smoothed to improve differentiability and enable gradient-based optimization.

- **Optimization workflows**:  
  `PVlib.jl` is intended to be embedded directly in Julia optimization workflows, including design optimization of photovoltaic systems and hybrid renewable energy platforms.

- **Floating PV utilities**:  
  `PVlib.jl` includes additional utilities for offshore and floating photovoltaic systems, including:
  - panel tilt and azimuth calculation from platform rotations,
  - rolling-average cell temperature for short-timescale simulations,
  - ocean surface albedo evaluation,
  - simplified geometric shading from nearby obstacles.

- **Verification against pvlib-python**:  
  Core photovoltaic modeling results from `PVlib.jl` have been compared against equivalent `pvlib-python` workflows. Small differences are expected because `PVlib.jl` smooths selected discontinuities to support differentiability.

In short, `pvlib-python` remains the more comprehensive photovoltaic modeling package, while `PVlib.jl` provides a Julia-native, differentiable, optimization-oriented implementation with additional support for floating PV applications.