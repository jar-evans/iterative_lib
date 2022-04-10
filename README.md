# iterative_lib

## Overview
Julia implementation of a bunch of iterative methods for solving linear systems. BIMs (basic iterative methods) aren't the most efficient, but provide a basis for more advanced approaches (preconditioners/smoothers)

### Current performance
![runtime](performance/runtime.png?raw=true "Runtime")

![memory_allocation](performance/memory_allocation.png?raw=true "Memory allocation")

## Problem statement
Solve for ![equation](https://latex.codecogs.com/svg.image?\textbf{u}):

![equation](https://latex.codecogs.com/svg.image?A\textbf{u}&space;=&space;\textbf{f})

At each iteration there exists the solution error (distance between current guess and exact solution):

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{e}^{k}=&space;\textbf{u}-\textbf{u}^{k})

However, if this is known then the solution is known, which it isn't. So the residual is defined:

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{r}^{k}=&space;\textbf{f}-A\textbf{u}^{k})

This is a measure of the aproximation of the current solution, and can then be used as a halting condition.

## Solution

Assuming that there exists a non-singular matrix ![equation](https://latex.codecogs.com/svg.image?\inline&space;M), a new matrix ![equation](https://latex.codecogs.com/svg.image?\inline&space;N) can be constructed:

![equation](https://latex.codecogs.com/svg.image?\inline&space;N=M-A)

This allows ![equation](https://latex.codecogs.com/svg.image?\inline&space;A\textbf{u}&space;=&space;\textbf{f}) to be rewritten to ![equation](https://latex.codecogs.com/svg.image?\inline&space;M\textbf{u}&space;=&space;N\textbf{u}+\textbf{f}), then multiplying through by ![equation](https://latex.codecogs.com/svg.image?\inline&space;M^{-1}) leads to an outline for an iterative scheme:

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{u}^{k&plus;1}&space;=&space;\textbf{u}^{k}&space;&plus;&space;M^{-1}\textbf{r}^k)

then reintroducing the definition of the residual leads to:

![equation](https://latex.codecogs.com/svg.image?\inline&space;\textbf{u}^{k&plus;1}&space;=&space;M^{-1}(N\textbf{u}^k&plus;\textbf{f}))

This equation acts as a basis for a range of stationary iterative solution methods, the choice of ![equation](https://latex.codecogs.com/svg.image?\inline&space;M) is what defines the method.

As iteration requires ... iteration, we want to make each update of the guess as cheap as possible, this is the case for diagonal or triangular matrices. With this in mind ![equation](https://latex.codecogs.com/svg.image?\inline&space;M) is commonly taken to be the diagonal or triangular components of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A). 

![equation](https://latex.codecogs.com/svg.image?\inline&space;A) is then interpreted as ![equation](https://latex.codecogs.com/svg.image?\inline&space;A&space;=&space;D&space;-&space;L&space;-&space;U)

![equation](https://latex.codecogs.com/svg.image?\inline&space;\xrightarrow[]{}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;D): diagonal of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A)

![equation](https://latex.codecogs.com/svg.image?\inline&space;\xrightarrow[]{}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;-L): lower triangle of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A)

![equation](https://latex.codecogs.com/svg.image?\inline&space;\xrightarrow[]{}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;-U): upper triangle of ![equation](https://latex.codecogs.com/svg.image?\inline&space;A)



### Jacobi
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;D) (![equation](https://latex.codecogs.com/svg.image?\inline&space;N&space;=&space;L+U))

 Which leads to:
 
![equation](https://latex.codecogs.com/svg.image?\textbf{u}^{k&plus;1}&space;=&space;D^{-1}[(L&plus;U)\textbf{u}^k&space;&plus;&space;\textbf{f}])
 
therefore producing the following update scheme:
 
![equation](https://latex.codecogs.com/svg.image?u^{k&plus;1}_{i}&space;=&space;\frac{1}{a_{ii}}[f_i&space;-&space;\sum_{j=1,j\neq&space;i}^{n}a_{ij}u_j^k]&space;&space;\forall_i&space;=&space;1,&space;...,&space;n)

As can be seen, the scheme is completely dependant on the previous iteration and independant of the other elements of the current iteration, making Jacobi natural to parallelise.
### Gauss-Seidal
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;D-L) (![equation](https://latex.codecogs.com/svg.image?\inline&space;N&space;=&space;U))

Which leads to:
 
![equation](https://latex.codecogs.com/svg.image?\textbf{u}^{k&plus;1}&space;=&space;(D&space;-&space;L)^{-1}(U\textbf{u}^k&space;&plus;&space;\textbf{f}))
 
therefore producing the following update scheme:

![equation](https://latex.codecogs.com/svg.image?u^{k&plus;1}_{i}&space;=&space;\frac{1}{a_{ii}}[f_i&space;-&space;\sum_{j=1}^{i-1}a_{ij}u_j^{k&plus;1}-&space;\sum_{j=i&plus;1}^{n}a_{ij}u_j^{k}]&space;&space;&space;&space;&space;\forall_i&space;=&space;1,&space;...,&space;n)

Unlike Jacobi, the scheme is dependant on the other elements of the current iteration, making Gauss-Seidal naturally serial, however ordering schemes such as red-black can be used to fix this.
### Damped Jacobi
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;\frac{1}{\omega}D)

For the damped case, the next guess becomes a weighted average of the current guess and the next guess as calculated by regular Jacobi:

![equation](https://latex.codecogs.com/svg.image?\textbf{u}^{k&plus;1}&space;=&space;(1-\omega)\textbf{u}^k&space;&plus;&space;\omega\textbf{u}^{k&plus;1}_{jacobi})
### Successive Over Relaxation (SOR)
![equation](https://latex.codecogs.com/svg.image?\inline&space;M&space;=&space;\frac{1}{\omega}D-L)

Similar to the damped Jacobi case, SOR is a weighted average:

![equation](https://latex.codecogs.com/svg.image?\textbf{u}^{k&plus;1}&space;=&space;(1-\omega)\textbf{u}^k&space;&plus;&space;\omega\textbf{u}^{k&plus;1}_{GS})
### Conjugate Gradient Descent

#TODO (actually a method that makes use of Krylov subspaces, but I thought it still fits here)

Sources:
<ul>
  <li> Richard S. Varga. Matrix Iterative Analysis. Springer, Berlin, New York, second edition, 2000 </li>
  <li> David M. Jr. Young. Iterative Solution of Large Linear Systems. Dover, Mineola, NY, USA, 2003 </li>
</ul>

