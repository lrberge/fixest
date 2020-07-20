
## Ongoing issues

#### Varying slopes convergence problem

There is a convergence problem that happens only when there are several variables with varying slopes (hereafter VVS). Having only one VVS is fine (note that I consider the slope fixed-effect as a variable too).

The issue is that I consider coefficients from variables with varying slopes as independent while they're in fact completely linked, and this *may* lead to convergence problems (note that convergence is usually fine though).

Now there's a warning popping when the convergence is not deemed good enough. 

The good news is that I realized there's a closed form solution to get all coefficients of VVS at once. But since it entails massive rewriting of the C++ code, I need to take at least 2/3 days of full time work to solve it and ensure everything is fine, and I don't know when that's going to happen. 

## Ongoing bugs



