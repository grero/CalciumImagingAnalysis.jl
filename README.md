# CalciumImagingAnalyses

Various tools to analyze calcium imaging traces in Julia.

## Installation

First install Julia by following these instructions: https://julialang.org/downloads/

Then, clone this repository, and start julia. Install all the dependencies by running

```julia
using Pkg
Pkg.activate(<path to repository>)
Pkg.instantiate()
```

where \<path to respository\> should be replaced with the actual path, for example "~/Documents/Julia/CalciumImagingAnalysis".

## Usage
See [this notebook](notebooks/usage.ipynb) for usage examples.

### From python
Thanks to the excellent PythonCall julia packace and the JuliaCall python package, you can call the functions from this Julia package also from python.

First, install the JuliaCall package into your current python envionment

```bash
pip install juliacall
```

Then, from python, first import the Main julia module, then import the Pkg module, and finally install this package

```python
from juliacall import Main as jl
jl.seval("using Pkg")
jl.Pkg.add(url="https://github.com/grero/CalciumImagingAnalysis.jl.git")
jl.seval("using CalciumImagingAnalyses") # import the module namespace
jl.seval("using CalciumImagingAnalyses: get_cell_data, get_behavioural_data, decode, train_decoder") # import useful functions
```

To carry on with the example from the notebook

```python
import os
import matplotlib.pylab as plt

datadir = "<path to data>"
cell_data_pre, timestamps_pre = jl.get_cell_data(os.path.join(datadir,"E10A","trace_10A_PRE3_unix.csv"),filter_accepted=False)
pokes_pre = get_behavioural_data(os.path.join(datadir, "E10A_E10B_Behavior/PRED3/E10A_PRED3.csv"))
q_pre = jl.decode(cell_data_pre, timestamps_pre, pokes_pre.rewarded_poke_time, pokes_pre.rewarded_poke_type,pratio=0.75)
plt.plot(q_pre)

```