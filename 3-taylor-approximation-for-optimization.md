# Lecture Notes: Taylor Approximation for Optimization

## Learning objectives

By the end of this lecture, students should be able to:

- Understand the geometric meaning of first- and second-order Taylor approximations.
- Derive the first-order approximation from the derivative.
- Derive the second-order approximation by studying how the slope changes.
- Understand why the factor $1/2$ appears in the second-order term.
- Extend the reasoning from one variable to two variables.
- Express Taylor approximations using the gradient and Hessian.
- Explain why Taylor approximation is a foundation for optimization methods.

# Part I. Taylor Approximation in One Variable

## 1. The basic question: How does a function change?

Consider a differentiable function

$$
f:\mathbb{R}\rightarrow\mathbb{R}.
$$

Suppose we know the value of the function at a point $a$, and we want to estimate its value at a nearby point $a+h$.

The change in the input is

$$
\Delta x=(a+h)-a=h.
$$

The corresponding change in the function is

$$
\Delta f=f(a+h)-f(a).
$$

Our objective is to understand how $\Delta f$ depends on the displacement $h$.

The fundamental idea is to use the information provided by the derivatives of $f$ at $a$.

## 2. First-order Taylor approximation

### 2.1 Starting from the derivative

The derivative at $a$ is defined as

$$
f'(a)=\lim_{h\rightarrow0}
\frac{f(a+h)-f(a)}{h}.
$$

For a sufficiently small displacement $h$, the ratio is approximately equal to $f'(a)$:

$$
\frac{f(a+h)-f(a)}{h}\approx f'(a).
$$

Multiplying both sides by $h$, we obtain

$$
f(a+h)-f(a)\approx f'(a)h.
$$

Therefore,

$$
\boxed{f(a+h)\approx f(a)+f'(a)h.}
$$

This is the first-order Taylor approximation.

### 2.2 Geometric interpretation

The derivative $f'(a)$ is the slope of the tangent line at $a$.

If the function is sufficiently smooth and $h$ is small, the tangent line provides a good local approximation to the function.

The tangent line is described by

$$
L(a+h)=f(a)+f'(a)h.
$$

Thus, the first-order approximation replaces the actual function with its tangent line near $a$.

The approximation has two components:

- $f(a)$: the initial function value.
- $f'(a)h$: the predicted change in the function.

The second component follows the principle

$$
\text{change in height}
\approx
\text{slope}\times\text{displacement}.
$$

### 2.3 Example

Consider

$$
f(x)=x^2
$$

at $a=2$.

We have

$$
f(2)=4,\qquad f'(x)=2x,\qquad f'(2)=4.
$$

For a displacement $h=0.1$, the first-order approximation gives

$$
f(2.1)\approx4+4(0.1)=4.4.
$$

The exact value is

$$
f(2.1)=2.1^2=4.41.
$$

The approximation is close because the displacement is small.

### 2.4 What does the first-order approximation neglect?

The tangent line uses only the slope at the starting point.

However, the actual slope may change as we move from $a$ to $a+h$.

Consequently, the function may gradually move away from the tangent line.

The first-order approximation neglects this variation in slope.

To improve the approximation, we must study how the slope itself changes.

# Part II. Second-Order Taylor Approximation in One Variable

## 3. From the slope to the change in slope

The first derivative $f'(x)$ describes the slope of the function.

The second derivative $f''(x)$ describes how that slope changes with position.

By definition,

$$
f''(a)=\lim_{h\rightarrow0}
\frac{f'(a+h)-f'(a)}{h}.
$$

For a small displacement $h$,

$$
\frac{f'(a+h)-f'(a)}{h}\approx f''(a).
$$

Multiplying by $h$,

$$
f'(a+h)-f'(a)\approx f''(a)h.
$$

Therefore,

$$
\boxed{f'(a+h)\approx f'(a)+f''(a)h.}
$$

This equation tells us that, when the second derivative is approximately constant over a small interval, the slope changes approximately linearly with the displacement.

The change in slope is

$$
\boxed{\Delta f'\approx f''(a)h.}
$$

This is the key to deriving the second-order approximation.

## 4. Deriving the second-order approximation using integration

### 4.1 Express the change in height as an integral

The exact change in the function is

$$
\Delta f=f(a+h)-f(a).
$$

By the fundamental theorem of calculus,

$$
\boxed{\Delta f=\int_a^{a+h}f'(x)\,dx.}
$$

This equation has a direct geometric interpretation: the change in the function equals the signed area under the derivative curve over the interval.

To calculate this area, we need to understand how the slope varies during the displacement.

### 4.2 Assume that the second derivative is constant

Suppose $f''(x)=c$ over the interval, where $c$ is constant.

Since

$$
\frac{df'}{dx}=f''(x)=c,
$$

integrating from $a$ to $a+s$ gives

$$
f'(a+s)-f'(a)=\int_0^s c\,dt.
$$

Therefore,

$$
f'(a+s)=f'(a)+cs.
$$

Since $c=f''(a)$,

$$
\boxed{f'(a+s)=f'(a)+f''(a)s.}
$$

The slope at displacement $s$ consists of the initial slope plus the change accumulated during that displacement.

### 4.3 Integrate the changing slope

Using $s$ as the displacement variable,

$$
\Delta f=\int_0^h f'(a+s)\,ds.
$$

Substitute the expression for the slope:

$$
\Delta f=\int_0^h
\left[f'(a)+f''(a)s\right]ds.
$$

Separate the integrals:

$$
\Delta f=
f'(a)\int_0^h ds+
f''(a)\int_0^h s\,ds.
$$

Evaluate the first integral:

$$
\int_0^h ds=h.
$$

Evaluate the second integral:

$$
\int_0^h s\,ds=\frac{h^2}{2}.
$$

Consequently,

$$
\boxed{\Delta f=f'(a)h+\frac12 f''(a)h^2.}
$$

Adding the initial function value gives

$$
\boxed{
f(a+h)=f(a)+f'(a)h+
\frac12 f''(a)h^2.
}
$$

Under the constant-second-derivative assumption, this expression is exact.

More generally, for a twice continuously differentiable function, the second-order Taylor approximation is

$$
\boxed{
f(a+h)\approx f(a)+f'(a)h+
\frac12 f''(a)h^2.
}
$$

The approximation error is $o(h^2)$ as $h\to0$.

## 5. Why does the factor $1/2$ appear?

The factor $1/2$ follows directly from integration.

The change in slope is approximately

$$
\Delta f'\approx f''(a)h.
$$

As we move through the interval, the slope correction grows approximately linearly from zero to $f''(a)h$.

The additional change in height is the integral of this growing correction:

$$
\Delta f_{\mathrm{correction}}
=\int_0^h f''(a)s\,ds.
$$

Therefore,

$$
\Delta f_{\mathrm{correction}}
=\frac12 f''(a)h^2.
$$

Geometrically, when $h>0$ and $f''(a)>0$, this correction is the area of a triangle whose base is $h$ and whose height is $f''(a)h$.

Thus,

$$
\text{additional height}
=\frac12\times
\text{base}\times\text{height}.
$$

The factor $1/2$ appears because the slope correction accumulates progressively rather than being equal to its final value throughout the movement.

## 6. Numerical example: comparing first- and second-order approximations

Consider

$$
f(x)=x^2
$$

at $a=2$, with $h=0.1$.

We have

$$
f(2)=4,\qquad f'(2)=4,\qquad f''(2)=2.
$$

The first-order approximation gives

$$
f(2.1)\approx4+4(0.1)=4.4.
$$

The second-order approximation gives

$$
\begin{aligned}
f(2.1)
&\approx4+4(0.1)+\frac12(2)(0.1)^2\\
&=4+0.4+0.01\\
&=4.41.
\end{aligned}
$$

The exact value is $4.41$.

In this example, the second-order approximation is exact because the function is quadratic and its second derivative is constant.

The difference between the two approximations is the curvature correction:

$$
\frac12 f''(2)h^2=0.01.
$$

This illustrates how the second derivative improves the prediction by accounting for the change in slope.

# Part III. Taylor Approximation in Two Variables

## 7. From one displacement to two displacements

Consider a function

$$
f:\mathbb{R}^2\rightarrow\mathbb{R}.
$$

We begin at a point

$$
P=(a,b)
$$

and move to

$$
Q=(a+h,b+k).
$$

The displacement is

$$
d=\begin{pmatrix}h\\k\end{pmatrix}.
$$

The exact change in the function is

$$
\boxed{\Delta f=f(a+h,b+k)-f(a,b).}
$$

There are now two coordinate displacements:

- $h$: displacement in the $x$-direction.
- $k$: displacement in the $y$-direction.

The objective is to understand how both displacements contribute to the change in the function.

## 8. Understanding the displacement geometrically

The total displacement can be decomposed into two coordinate movements:

$$
\begin{pmatrix}h\\k\end{pmatrix}
=
\begin{pmatrix}h\\0\end{pmatrix}
+
\begin{pmatrix}0\\k\end{pmatrix}.
$$

We can therefore reach $Q$ by:

1. Moving horizontally by $h$, from $(a,b)$ to $(a+h,b)$.
2. Moving vertically by $k$, from $(a+h,b)$ to $(a+h,b+k)$.

Alternatively, we can move directly along the vector $(h,k)$.

The final position is the same because vector displacements add.

However, the height changes during these movements depend on the slopes of the surface, and these slopes can vary with position.

This distinction between displacement and height change is fundamental:

- The displacement in the domain is a vector.
- The change in the function is a scalar.
- The partial derivatives determine how the displacement affects the height.

## 9. First-order approximation in two variables

### 9.1 Partial derivatives as slopes

The partial derivative

$$
f_x(a,b)
$$

measures the slope of the surface in the $x$-direction while keeping $y=b$ fixed.

Similarly,

$$
f_y(a,b)
$$

measures the slope in the $y$-direction while keeping $x=a$ fixed.

For small coordinate displacements, the corresponding first-order height contributions are

$$
\Delta f_x\approx f_x(a,b)h
$$

and

$$
\Delta f_y\approx f_y(a,b)k.
$$

Adding them gives

$$
\boxed{
\Delta f\approx f_x(a,b)h+f_y(a,b)k.
}
$$

Therefore,

$$
\boxed{
f(a+h,b+k)\approx
f(a,b)+f_x(a,b)h+f_y(a,b)k.
}
$$

This is the first-order Taylor approximation in two variables.

### 9.2 Geometric interpretation

The first-order approximation uses the two initial slopes to predict the height change.

It is the local tangent-plane approximation to the surface.

The horizontal displacement contributes approximately $f_xh$, and the vertical displacement contributes approximately $f_yk$.

The sum gives the first-order prediction of the total change in height.

However, just as in one variable, this approximation assumes that the initial slopes remain representative during the movement.

To improve the approximation, we must study how both slopes change.

# Part IV. Second-Order Taylor Approximation in Two Variables

## 10. How does the horizontal slope change?

The horizontal slope is

$$
f_x(x,y).
$$

Its value depends on both coordinates. Therefore, moving horizontally or vertically can change it.

There are two contributions.

### 10.1 Change caused by horizontal displacement

The rate at which $f_x$ changes when $x$ changes is

$$
f_{xx}(a,b).
$$

For a small horizontal displacement $h$,

$$
\boxed{
\Delta f_x\big|_{\mathrm{horizontal}}
\approx f_{xx}(a,b)h.
}
$$

The second partial derivative $f_{xx}$ measures the change in the horizontal slope as we move horizontally.

### 10.2 Change caused by vertical displacement

The rate at which $f_x$ changes when $y$ changes is

$$
f_{xy}(a,b).
$$

For a small vertical displacement $k$,

$$
\boxed{
\Delta f_x\big|_{\mathrm{vertical}}
\approx f_{xy}(a,b)k.
}
$$

The mixed partial derivative $f_{xy}$ measures how the horizontal slope changes when we move vertically.

### 10.3 Total change in the horizontal slope

Adding the two contributions gives

$$
\boxed{
\Delta f_x\approx
f_{xx}(a,b)h+f_{xy}(a,b)k.
}
$$

The horizontal slope changes because of both coordinate displacements.

## 11. How does the vertical slope change?

The vertical slope is

$$
f_y(x,y).
$$

Again, it can change because of movement in either coordinate.

### 11.1 Change caused by horizontal displacement

The rate at which $f_y$ changes horizontally is

$$
f_{yx}(a,b).
$$

Therefore,

$$
\boxed{
\Delta f_y\big|_{\mathrm{horizontal}}
\approx f_{yx}(a,b)h.
}
$$

### 11.2 Change caused by vertical displacement

The rate at which $f_y$ changes vertically is

$$
f_{yy}(a,b).
$$

Therefore,

$$
\boxed{
\Delta f_y\big|_{\mathrm{vertical}}
\approx f_{yy}(a,b)k.
}
$$

### 11.3 Total change in the vertical slope

Adding the contributions,

$$
\boxed{
\Delta f_y\approx
f_{yx}(a,b)h+f_{yy}(a,b)k.
}
$$

We now have expressions for how both initial slopes change during the displacement.

## 12. How do the changes in slope affect the height?

Recall the first-order approximation:

$$
\Delta f\approx f_xh+f_yk.
$$

The slopes in this expression are evaluated at the starting point $P=(a,b)$.

However, as we move toward $Q$, the slopes change. We must therefore account for the additional height change caused by their variation.

To derive the second-order formula rigorously, we integrate the changing slopes along the straight-line path from $P$ to $Q$.

### 12.1 Parameterize the displacement

Define

$$
x(t)=a+th,\qquad y(t)=b+tk,
\qquad 0\leq t\leq1.
$$

At $t=0$, we are at $P$. At $t=1$, we are at $Q$.

Along this path,

$$
\frac{dx}{dt}=h,\qquad
\frac{dy}{dt}=k.
$$

By the chain rule,

$$
\frac{df}{dt}
=f_x(x(t),y(t))h+
f_y(x(t),y(t))k.
$$

Integrating from $t=0$ to $t=1$,

$$
\boxed{
\Delta f=\int_0^1
\left[f_x(x(t),y(t))h+
f_y(x(t),y(t))k\right]dt.
}
$$

This expression gives the exact change in height along the path.

### 12.2 Approximate the changing slopes

The total change in $f_x$ along the path is approximately

$$
\Delta f_x\approx f_{xx}h+f_{xy}k.
$$

Similarly,

$$
\Delta f_y\approx f_{yx}h+f_{yy}k.
$$

Assume that the second derivatives are constant along the path. The slopes then vary linearly with $t$:

$$
\begin{aligned}
f_x(t)&=f_x(P)+t(f_{xx}h+f_{xy}k),\\
f_y(t)&=f_y(P)+t(f_{yx}h+f_{yy}k).
\end{aligned}
$$

These equations express the changing slopes in terms of their initial values and the accumulated changes.

### 12.3 Substitute into the height integral

Substituting into the integral gives

$$
\begin{aligned}
\Delta f=\int_0^1 \Big[&
f_x(P)h+f_y(P)k\\
&+t(f_{xx}h+f_{xy}k)h\\
&+t(f_{yx}h+f_{yy}k)k
\Big]dt.
\end{aligned}
$$

Separate the initial contributions from the slope-change contributions:

$$
\begin{aligned}
\Delta f={}&
f_x(P)h+f_y(P)k\\
&+(f_{xx}h^2+f_{xy}hk)\int_0^1t\,dt\\
&+(f_{yx}hk+f_{yy}k^2)\int_0^1t\,dt.
\end{aligned}
$$

Since

$$
\int_0^1t\,dt=\frac12,
$$

we obtain

$$
\begin{aligned}
\Delta f={}&f_x(P)h+f_y(P)k\\
&+\frac12 f_{xx}(P)h^2\\
&+\frac12 f_{xy}(P)hk\\
&+\frac12 f_{yx}(P)hk\\
&+\frac12 f_{yy}(P)k^2.
\end{aligned}
$$

This is the second-order expansion before combining the mixed terms.

## 13. Combine the mixed terms

If the second partial derivatives are continuous near $P$, then Clairaut's theorem gives

$$
f_{xy}(P)=f_{yx}(P).
$$

Therefore,

$$
\frac12 f_{xy}hk+\frac12 f_{yx}hk
=f_{xy}hk.
$$

The final expanded formula is

$$
\boxed{
\begin{aligned}
f(a+h,b+k)\approx{}&f(a,b)\\
&+f_x(a,b)h+f_y(a,b)k\\
&+\frac12 f_{xx}(a,b)h^2\\
&+f_{xy}(a,b)hk\\
&+\frac12 f_{yy}(a,b)k^2.
\end{aligned}
}
$$

Equivalently, the change in the function is approximated by

$$
\boxed{
\begin{aligned}
\Delta f\approx{}&
f_xh+f_yk\\
&+\frac12 f_{xx}h^2
+f_{xy}hk
+\frac12 f_{yy}k^2.
\end{aligned}
}
$$

All derivatives are evaluated at $(a,b)$.

The approximation error is $o(h^2+k^2)$ as $(h,k)\to(0,0)$, assuming $f$ is twice continuously differentiable near the starting point.

If the second derivatives are constant, the expression is exact along the straight-line path.

## 14. Understanding the terms geometrically

The expanded formula contains three types of contributions.

### Initial slopes

$$
f_xh+f_yk.
$$

These describe the height change predicted using the slopes at the starting point.

### Curvature along the coordinate directions

$$
\frac12 f_{xx}h^2
\qquad\text{and}\qquad
\frac12 f_{yy}k^2.
$$

These account for the changes in the horizontal and vertical slopes caused by movement along their respective coordinate directions.

### Interaction between the two displacements

$$
f_{xy}hk.
$$

This term accounts for how the two coordinate directions interact through the changing slopes.

If either $h=0$ or $k=0$, the mixed contribution vanishes.

The complete formula combines the initial height prediction, the curvature corrections in each coordinate direction, and the interaction between the two directions.

# Part V. General Form Using the Gradient and Hessian

## 15. From two coordinates to a vector

The same reasoning applies to a function with $n$ input variables:

$$
f:\mathbb{R}^n\rightarrow\mathbb{R}.
$$

Let the current point be

$$
x=
\begin{pmatrix}
x_1\\
x_2\\
\vdots\\
x_n
\end{pmatrix}
$$

and let the displacement be

$$
d=
\begin{pmatrix}
d_1\\
d_2\\
\vdots\\
d_n
\end{pmatrix}.
$$

The new point is

$$
x+d.
$$

We seek an approximation to $f(x+d)$ using the derivatives evaluated at $x$.

## 16. The gradient: collecting the first derivatives

The gradient is the vector of first partial derivatives:

$$
\boxed{
\nabla f(x)=
\begin{pmatrix}
\frac{\partial f}{\partial x_1}\\
\frac{\partial f}{\partial x_2}\\
\vdots\\
\frac{\partial f}{\partial x_n}
\end{pmatrix}.
}
$$

The first-order contribution is the dot product of the gradient and the displacement:

$$
\boxed{\nabla f(x)^Td.}
$$

Expanding the dot product,

$$
\nabla f(x)^Td
=\sum_{i=1}^{n}
\frac{\partial f}{\partial x_i}d_i.
$$

Each partial derivative gives the rate of change along one coordinate direction, and multiplication by the corresponding displacement gives its first-order contribution.

Adding the contributions gives the first-order approximation:

$$
\boxed{
f(x+d)\approx f(x)+\nabla f(x)^Td.
}
$$

## 17. The Hessian: collecting the second derivatives

The Hessian is the matrix of second partial derivatives:

$$
\boxed{
H_f(x)=
\begin{pmatrix}
\frac{\partial^2 f}{\partial x_1^2}
&
\frac{\partial^2 f}{\partial x_1\partial x_2}
&\cdots&
\frac{\partial^2 f}{\partial x_1\partial x_n}
\\
\frac{\partial^2 f}{\partial x_2\partial x_1}
&
\frac{\partial^2 f}{\partial x_2^2}
&\cdots&
\frac{\partial^2 f}{\partial x_2\partial x_n}
\\
\vdots&\vdots&\ddots&\vdots\\
\frac{\partial^2 f}{\partial x_n\partial x_1}
&
\frac{\partial^2 f}{\partial x_n\partial x_2}
&\cdots&
\frac{\partial^2 f}{\partial x_n^2}
\end{pmatrix}.
}
$$

Each entry describes how a component of the gradient changes when we move along a coordinate direction.

Multiplying the Hessian by the displacement gives the first-order change in the gradient:

$$
\boxed{\Delta(\nabla f)\approx H_f(x)d.}
$$

Taking the dot product with the displacement once more gives

$$
d^TH_f(x)d.
$$

This is a scalar that describes the second-order contribution along the chosen displacement.

The factor $1/2$ arises from integrating the linearly accumulating change in the gradient along the path.

Therefore, the second-order Taylor approximation in $n$ variables is

$$
\boxed{
f(x+d)\approx
f(x)+\nabla f(x)^Td+
\frac12 d^TH_f(x)d.
}
$$

This is the general matrix form of the expanded formula derived for two variables.

## 18. Final connection to optimization

Taylor approximation provides a local model of an objective function near a current point $x$.

The first-order model is

$$
m_1(d)=f(x)+\nabla f(x)^Td.
$$

It uses the gradient to predict the change in the objective.

The second-order model is

$$
m_2(d)=f(x)+\nabla f(x)^Td+
\frac12d^TH_f(x)d.
$$

It also accounts for changes in slope through the Hessian.

These models lead naturally to optimization methods:

- **Gradient descent** uses the gradient to choose a descent direction.
- **Line search** chooses a step length along a selected direction.
- **Newton's method** uses the quadratic model and, when the Hessian is invertible, obtains a step by solving
  $$
  H_f(x)d=-\nabla f(x).
  $$

The central principle is that optimization methods use local information about the objective function to choose useful displacements.

Taylor approximation explains how that local information is constructed: first from the slopes, then from the changes in slopes, and finally from the gradient and Hessian.
