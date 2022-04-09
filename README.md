# iterative_lib

## Overview
Julia implementation of a bunch of iterative methods for solving linear systems. BIMs (basic iterative methods) aren't the most efficient, but provide a basis for more advanced approaches (preconditioners/smoothers)

![runtime](performance/runtime.png?raw=true "Runtime")

![memory_allocation](performance/memory_allocation.png?raw=true "Memory allocation")
## Problem statement
Solve:
![equation](https://latex.codecogs.com/svg.image?A\textbf{u}&space;=&space;\textbf{f})

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{e}^{k}=&space;\textbf{u}-\textbf{u}^{k})

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{r}^{k}=&space;\textbf{f}-A\textbf{u}^{k})

Assuming that there exists a non-singular matrix ![equation](https://latex.codecogs.com/svg.image?\inline&space;M), a new matrix ![equation](https://latex.codecogs.com/svg.image?\inline&space;N) can be constructed:

![equation](https://latex.codecogs.com/svg.image?\inline&space;N=M-A)

This allows ![equation](https://latex.codecogs.com/svg.image?\inline&space;A\textbf{u}&space;=&space;\textbf{f}) to be rewritten to ![equation](https://latex.codecogs.com/svg.image?\inline&space;M\textbf{u}&space;=&space;N\textbf{u}+\textbf{f}), then multiplying through by ![equation](https://latex.codecogs.com/svg.image?\inline&space;M^{-1}) leads to an outline for an iterative scheme:

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{u}^{k&plus;1}&space;=&space;\textbf{u}^{k}&space;&plus;&space;M^{-1}\textbf{r}^k)

then reintroducing the definition of the residual leads to:

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{u}^{k&plus;1}&space;=&space;M^{-1}(N\textbf{u}^k&plus;\textbf{f}))

This equation acts as a basis for a range of stationary iterative solution methods, the choice of ![equation](https://latex.codecogs.com/svg.image?\inline&space;M) is what defines the method.

As iteration requires ... iteration, we want to make each update of the guess as cheap as possible, this is the case for diagonal or triangular matrices. With this in mind ![equation](https://latex.codecogs.com/svg.image?\inline&space;M) is commonly taken to be the diagonal or triangular components of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A). 

![equation](https://latex.codecogs.com/svg.image?\inline&space;A) is then interpreted as ![equation](https://latex.codecogs.com/svg.image?\inline&space;A&space;=&space;D&space;-&space;E&space;-&space;F&space;\in&space;\mathbb{R}^{n\times&space;n})

![equation](https://latex.codecogs.com/svg.image?\inline&space;\xrightarrow[]{}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;D): diagonal of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A)

![equation](https://latex.codecogs.com/svg.image?\inline&space;\xrightarrow[]{}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;-E): lower triangle of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A)

![equation](https://latex.codecogs.com/svg.image?\inline&space;\xrightarrow[]{}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;-F): upper triangle of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A)



### Jacobi
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;D) (![equation](https://latex.codecogs.com/svg.image?\inline&space;N&space;=&space;E+F))

### Gauss-Seidal
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;D-E) (![equation](https://latex.codecogs.com/svg.image?\inline&space;N&space;=&space;F))
### Damped Jacobi
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;\frac{1}{\omega}D)
### Successive Over Relaxation (SOR)
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;\frac{1}{\omega}D-E)
### Conjugate Gradient Descent

Sources:
<ul>
  <li> Richard S. Varga. Matrix Iterative Analysis. Springer, Berlin, New York, second edition, 2000 </li>
  <li> David M. Jr. Young. Iterative Solution of Large Linear Systems. Dover, Mineola, NY, USA, 2003 </li>
</ul>

