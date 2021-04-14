# VQTPlot

Generates the plots of the workload distributions of a special kind of queueing system.

This queueing system has `c` servers and an arrival stream of customers of Poisson type with homogeneous parameter `λ`.

Whenever a customer arrives it is tagged as customer of type 1 on arrival the workload is below the threshold `k`, otherwise it is tagged as customer of type 2.

All customers have a service time exponential distributed, type 1 customers with parameter `μ1` and type 2 customers with parameter `μ2`.


You can run the code in your browser here on [GitHub Pages](https://brdauria.github.io/workloadPlot/)

or run the jupyther notebook in python on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brdauria/workloadPlot/HEAD)
