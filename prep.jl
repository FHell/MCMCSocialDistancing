# This file instantiates the committed Manifest.toml

import Pkg
Pkg.activate(".")
Pkg.instantiate()

using LightGraphs
using Statistics
using BSON
