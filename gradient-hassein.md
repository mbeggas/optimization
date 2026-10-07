# Understanding $v^T H v$: Directional Curvature and the Hessian

## 1. From a point and a direction to a line

Let

$$p \in \mathbb{R}^n$$

be a point, represented by its position vector, and let

$$v \in \mathbb{R}^n$$

be a direction vector.

Consider the expression

$$\boxed{x(t)=p+tv,\qquad t\in\mathbb{R}.}$$

This describes a straight line passing through $p$ in the direction $v$.

### Why?

For two values $t_1$ and $t_2$,

$$x(t_1)=p+t_1v,$$

and

$$x(t_2)=p+t_2v.$$

Their difference is

$$x(t_2)-x(t_1) = (p+t_2v)-(p+t_1v) = (t_2-t_1)v.$$

Thus the displacement between any two points on this trajectory is a scalar multiple of $v$.

Therefore,

$$\boxed{\text{the line has direction }v.}$$

The vector $p$ determines **where the line is located**, while $v$ determines **its direction**.

---

## 2. Restricting a multivariable function to a line

Suppose

$$f:\mathbb{R}^n\rightarrow\mathbb{R}.$$

The function $f$ describes a potentially high-dimensional landscape.

Instead of studying the whole landscape at once, suppose we start at $p$ and move only along the direction $v$:

$$x(t)=p+tv.$$

We can then observe the value of $f$ along this line:

$$\boxed{g(t)=f(p+tv).}$$

This is a very important transformation.

Originally,

$$f(x)$$

is a function of $n$ variables.

After restricting ourselves to the line $p+tv$,

$$g(t)=f(p+tv)$$

is a function of only **one variable**, $t$.

Thus,

$$\boxed{\text{multivariable function} \quad\longrightarrow\quad \text{one-dimensional function along a chosen direction}}$$

This allows us to use ordinary one-dimensional derivatives to study the geometry of $f$ in a particular direction.

---

## 3. The first derivative: slope in a direction

We have

$$g(t)=f(p+tv).$$

Differentiate with respect to $t$.

By the chain rule,

$$g'(t) = \nabla f(p+tv)^T v.$$

At $t=0$,

$$\boxed{g'(0)=\nabla f(p)^T v.}$$

This quantity has a simple geometric interpretation:

$$\boxed{\nabla f(p)^T v = \text{rate of change of }f\text{ at }p\text{ in direction }v.}$$

Thus, the gradient contains information about the slope of the function in every possible direction.

For a particular direction $v$, the scalar

$$\nabla f(p)^Tv$$

extracts the slope in that direction.

---

## 4. The second derivative: curvature in a direction

The first derivative tells us how quickly the function is changing.

The second derivative tells us how this rate of change itself changes.

Starting from

$$g'(t)=\nabla f(p+tv)^Tv,$$

differentiate again:

$$g''(t) = \frac{d}{dt} \left[ \nabla f(p+tv)^Tv \right].$$

The derivative of the gradient is the Hessian:

$$H(x)=\nabla^2 f(x).$$

Therefore,

$$\boxed{g''(t)=v^T H(p+tv)v.}$$

In particular, at $t=0$,

$$\boxed{g''(0)=v^T H(p)v.}$$

This gives the fundamental interpretation:

> **$v^T H(p)v$ is the second derivative of $f$ at $p$ when we move in the direction $v$.**

Equivalently,

$$\boxed{v^T H(p)v = \text{directional curvature of }f\text{ at }p\text{ along }v.}$$

This is the central meaning of the quadratic form $v^THv$.

---

## 5. What exactly happens in $v^T H v$?

Let

$$H = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix}$$

and

$$v = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}.$$

The expression

$$v^THv$$

contains two matrix-vector operations.

### Step 1: Apply the Hessian to $v$

First calculate

$$Hv.$$

We obtain

$$Hv = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = \begin{pmatrix} 2v_1 \\ 4v_2 \end{pmatrix}.$$

The result is still a vector.

### Step 2: Project onto $v$

Now multiply by $v^T$:

$$v^THv = \begin{pmatrix} v_1 & v_2 \end{pmatrix} \begin{pmatrix} 2v_1 \\ 4v_2 \end{pmatrix}.$$

Therefore,

$$\boxed{v^THv=2v_1^2+4v_2^2.}$$

The result is a scalar.

That scalar is precisely the curvature of the function along the direction $v$.

---

## 6. Why do we use $v$ twice?

This is an important conceptual question.

Why do we calculate

$$v^THv$$

rather than simply

$$Hv?$$

The two expressions have different meanings.

### $Hv$

is a vector.

It describes how the Hessian acts on the direction $v$.

### $v^THv$

is a scalar.

It measures the second-order effect associated with moving along $v$.

Conceptually,

$$\boxed{v \overset{H}{\longrightarrow} Hv \overset{v^T}{\longrightarrow} v^THv.}$$

The Hessian describes the second-order behavior of the function, and the quadratic form $v^THv$ extracts the part of that behavior associated with the chosen direction $v$.

Most importantly, this is not merely an interpretation imposed on the formula:

$$\boxed{v^THv=g''(0)}$$

follows directly from the chain rule.

---

## 7. Example: $f(x,y)=x^2+2y^2$

Consider

$$f(x,y)=x^2+2y^2.$$

Its gradient is

$$\nabla f(x,y) = \begin{pmatrix} 2x \\ 4y \end{pmatrix},$$

and its Hessian is

$$\boxed{H = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix}.}$$

We can now examine the curvature in different directions.

---

### 7.1 Moving horizontally

Choose

$$v = \begin{pmatrix} 1 \\ 0 \end{pmatrix}.$$

This means that we move in the $x$-direction.

The directional curvature is

$$v^THv = \begin{pmatrix} 1 & 0 \end{pmatrix} \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix} \begin{pmatrix} 1 \\ 0 \end{pmatrix} = 2.$$

We can verify this directly.

Along the $x$-direction,

$$y=0,$$

so

$$f(x,0)=x^2.$$

Its second derivative is

$$\frac{d^2}{dx^2}x^2=2.$$

Thus,

$$\boxed{v^THv=2.}$$

---

### 7.2 Moving vertically

Now choose

$$v = \begin{pmatrix} 0 \\ 1 \end{pmatrix}.$$

We move in the $y$-direction.

Then

$$v^THv=4.$$

Again, verify this directly.

Along the $y$-direction,

$$x=0,$$

so

$$f(0,y)=2y^2.$$

Therefore,

$$\frac{d^2}{dy^2}(2y^2)=4.$$

Hence,

$$\boxed{v^THv=4.}$$

The surface has different curvature in different directions.

---

## 8. Moving diagonally

Now choose

$$v = \begin{pmatrix} 1 \\ 1 \end{pmatrix}.$$

This means that we move along the diagonal direction $y=x$.

The quadratic form gives

$$v^THv = 2(1)^2+4(1)^2 = 6.$$

Let's verify this directly.

Starting from the origin,

$$p = \begin{pmatrix} 0 \\ 0 \end{pmatrix},$$

we move along

$$x(t)=p+tv = \begin{pmatrix} t \\ t \end{pmatrix}.$$

Therefore,

$$f(t,t) = t^2+2t^2 = 3t^2.$$

Its second derivative is

$$\frac{d^2}{dt^2}(3t^2)=6.$$

Thus,

$$\boxed{v^THv=6.}$$

This confirms that $v^THv$ is not an arbitrary matrix calculation.

It reproduces the ordinary one-dimensional second derivative after the multidimensional function has been restricted to a line.

---

## 9. Why $v^THv$ can be called curvature

In one dimension, the second derivative

$$f''(x)$$

measures the local curvature of a function.

In several dimensions, there is no single curvature value at a point.

Why?

Because the surface can curve differently in different directions.

Therefore, at a point $p$, we need to ask:

> What is the curvature if I move in this particular direction $v$?

The answer is

$$\boxed{v^TH(p)v.}$$

Thus, instead of thinking of the Hessian as providing one curvature number, think of it as providing a **curvature value for every possible direction**.

$$\boxed{v \quad\longrightarrow\quad v^TH(p)v}$$

means:

$$\boxed{\text{direction} \quad\longrightarrow\quad \text{curvature in that direction}.}$$

---

## 10. Positive semidefinite Hessian

We can now understand the definition of positive semidefiniteness geometrically.

A symmetric matrix $H$ is positive semidefinite if

$$\boxed{v^THv\geq0 \qquad\forall v.}$$

But we have just established that

$$v^THv$$

is the curvature in direction $v$.

Therefore,

$$\boxed{H\succeq0 \quad\Longleftrightarrow\quad \text{the curvature is nonnegative in every direction}.}$$

In geometric language:

> **The function does not bend downward in any direction.**

Some directions may have positive curvature, while others may have zero curvature.

---

## 11. Positive definite Hessian

If instead

$$\boxed{v^THv>0 \qquad\forall v\neq0,}$$

then the Hessian is positive definite:

$$\boxed{H\succ0.}$$

Geometrically:

> **The function has strictly positive curvature in every nonzero direction.**

There are no flat directions.

For example,

$$H = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix}$$

is positive definite because

$$v^THv = 2v_1^2+4v_2^2>0$$

for every

$$v\neq0.$$

---

## 12. What happens when $v^THv<0$?

Consider

$$f(x,y)=x^2-y^2.$$

Its Hessian is

$$\boxed{H = \begin{pmatrix} 2 & 0 \\ 0 & -2 \end{pmatrix}.}$$

Choose

$$v = \begin{pmatrix} 0 \\ 1 \end{pmatrix}.$$

Then

$$v^THv=-2.$$

Therefore, the curvature in the $y$-direction is negative.

Indeed,

$$f(0,y)=-y^2,$$

which curves downward.

But in the $x$-direction,

$$f(x,0)=x^2,$$

which curves upward.

Therefore, the function has:

* positive curvature in some directions;
* negative curvature in other directions.

This is the geometric signature of a **saddle**.

Algebraically,

$$\boxed{v^THv = 2v_1^2-2v_2^2,}$$

which can be either positive or negative depending on $v$.

Thus the Hessian is **indefinite**.

---

## 13. The complete geometric chain

We can now connect all the concepts.

Start with a point $p$ and a direction $v$:

$$\boxed{x(t)=p+tv.}$$

This gives a straight line through $p$.

Restrict the function to this line:

$$\boxed{g(t)=f(p+tv).}$$

The first derivative is

$$\boxed{g'(0)=\nabla f(p)^Tv,}$$

which is the **slope in direction $v$**.

The second derivative is

$$\boxed{g''(0)=v^TH(p)v,}$$

which is the **curvature in direction $v$**.

Therefore,

$$\boxed{v^TH(p)v\geq0}$$

means:

> The function curves upward when leaving $p$ in direction $v$.

If this holds for every direction,

$$\boxed{v^TH(p)v\geq0 \qquad\forall v,}$$

then:

> The function has nonnegative curvature in every possible direction.

This is exactly the geometric meaning of

$$\boxed{\nabla^2f(p)\succeq0.}$$

---

## 14. From local curvature to convexity

If the Hessian is positive semidefinite at every point,

$$\boxed{\nabla^2f(x)\succeq0 \qquad\forall x,}$$

then the function is convex (under the usual differentiability assumptions on a convex domain).

Geometrically:

$$\boxed{\text{PSD Hessian} \Rightarrow \text{nonnegative curvature in every direction} \Rightarrow \text{convex function}.}$$

The important conceptual point is that **convexity is not fundamentally about matrices**.

The matrix is a compact way of encoding a geometric fact:

> **At every point, and in every direction, the function does not curve downward.**

The quadratic form

$$\boxed{v^THv}$$

is the mechanism that extracts the curvature corresponding to one particular direction.
