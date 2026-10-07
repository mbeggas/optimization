# From Taylor Approximation to Gradient Descent and Line Search

## 1. The basic problem: how can we minimize a function?

Suppose we want to solve

$$
\min_{x\in\mathbb{R}^n} f(x).
$$

We are at the current point

$$
x_k.
$$

The question is:

> **Where should we move from $x_k$ to make $f$ smaller?**

We could try to understand the entire function $f$, but optimization methods usually work with **local information** around the current point.

This is where Taylor approximation enters.

---

# 2. Taylor approximation: replacing a complicated function by a local approximation

## 2.1 One-dimensional case

For a function $f(x)$, Taylor expansion around $x_0$ gives

$$
f(x_0+h)
\approx
f(x_0)
+
f'(x_0)h
+
\frac12 f''(x_0)h^2.
$$

There are three levels of approximation.

### Zeroth order

$$
f(x_0+h)\approx f(x_0)
$$

We know only the current value.

### First order

$$
\boxed{
f(x_0+h)
\approx
f(x_0)+f'(x_0)h
}
$$

We know the value and the slope.

This is the **tangent-line approximation**.

### Second order

$$
\boxed{
f(x_0+h)
\approx
f(x_0)+f'(x_0)h+\frac12f''(x_0)h^2
}
$$

Now we also know the local curvature.

---

# 3. Taylor approximation in several dimensions

For

$$
f:\mathbb{R}^n\rightarrow\mathbb{R},
$$

the displacement from $x$ is a vector $s$.

The first-order Taylor approximation is

$$
\boxed{
f(x+s)
\approx
f(x)+\nabla f(x)^Ts
}
$$

where

$$
\nabla f(x)=
\begin{pmatrix}
\frac{\partial f}{\partial x_1}\\
\vdots\\
\frac{\partial f}{\partial x_n}
\end{pmatrix}.
$$

The second-order approximation is

$$
\boxed{
f(x+s)
\approx
f(x)
+
\nabla f(x)^Ts +
\frac12s^T H(x)s
}
$$

where $H(x)$ is the Hessian.

---

# 4. What does $\nabla f(x)^Ts$ mean?

This is important for understanding optimization.

Suppose we choose a direction $v$ and move along it:

$$
x(t)=x+tv.
$$

We can look at the function only along this line:

$$
g(t)=f(x+tv).
$$

Then

$$
g'(0)=\nabla f(x)^Tv.
$$

Therefore,

$$
\boxed{\nabla f(x)^Tv}
$$

is the **slope of the function in direction $v$**.

So Taylor becomes

$$
f(x+tv)
\approx
f(x)
+
t\nabla f(x)^Tv.
$$

This gives us a very simple question:

> **Which direction should we choose to make this quantity as negative as possible?**

---

# 5. Finding a direction that decreases the function

We want

$$
\nabla f(x)^Tv<0.
$$

For example, if

$$
v=-\nabla f(x),
$$

then

$$
\nabla f(x)^Tv
=
\nabla f(x)^T[-\nabla f(x)]
$$

and therefore

$$
\boxed{
\nabla f(x)^Tv =
-\|\nabla f(x)\|^2<0
}
$$

whenever

$$
\nabla f(x)\neq0.
$$

Therefore,

$$
f(x+tv)<f(x)
$$

for sufficiently small positive $t$.

This gives us the fundamental idea of **gradient descent**.

---

# 6. Why is the negative gradient special?

We can ask an even stronger question.

Suppose we restrict ourselves to unit directions:

$$
\|v\|=1.
$$

We want the direction that produces the greatest instantaneous decrease:

$$
\min_{\|v\|=1}\nabla f(x)^Tv.
$$

By the Cauchy–Schwarz inequality,

$$
\nabla f(x)^Tv
\geq
-\|\nabla f(x)\|\|v\|.
$$

Since $\|v\|=1$,

$$
\nabla f(x)^Tv
\geq
-\|\nabla f(x)\|.
$$

Equality occurs when

$$
\boxed{
v=-\frac{\nabla f(x)}{\|\nabla f(x)\|}
}
$$

Therefore:

$$
\boxed{
-\nabla f(x)
\text{ is the direction of steepest descent.}
}
$$

This is the geometric origin of gradient descent.

---

# 7. From direction to an actual optimization step

Knowing the direction is not enough.

Suppose we know that

$$
d_k=-\nabla f(x_k)
$$

is a good direction.

We still need to decide:

> **How far should we move?**

Let the step length be $\alpha_k>0$.

Then

$$
\boxed{
x_{k+1}=
x_k+\alpha_kd_k
}
$$

and since

$$
d_k=-\nabla f(x_k),
$$

we obtain

$$
\boxed{
x_{k+1}
=
x_k-\alpha_k\nabla f(x_k)
}
$$

This is the **gradient descent update**.

Notice that two decisions are involved:

$$
\boxed{
\text{direction}
\quad+\quad
\text{step length}
}
$$

Gradient descent determines the direction from the gradient.

We still need a method for determining the step length.

This leads to **line search**.

---

# 8. Line search: reduce the multidimensional problem to one dimension

Suppose we are at $x_k$ and have chosen a direction $d_k$.

Consider all points obtained by moving along this direction:

$$
\boxed{
x_k+\alpha d_k
}
$$

where $\alpha$ is a scalar.

This is a line in the $n$-dimensional space.

Now define a new one-dimensional function:

$$
\boxed{
\phi(\alpha)=
f(x_k+\alpha d_k)
}
$$

Instead of asking

$$
\min_x f(x),
$$

we temporarily ask

$$
\boxed{
\min_{\alpha\geq0}\phi(\alpha).
}
$$

This is the fundamental idea of **line search**.

We have transformed:

$$
\text{multidimensional optimization}
$$

into

$$
\text{one-dimensional optimization along a chosen direction}.
$$

---

# 9. How does Taylor approximation describe the line search?

Apply Taylor expansion to

$$
\phi(\alpha)=
f(x_k+\alpha d_k).
$$

Around $\alpha=0$,

$$
\phi(\alpha)
\approx
\phi(0)
+
\alpha\phi'(0)
+
\frac12\alpha^2\phi''(0).
$$

Since

$$
\phi(0)=f(x_k),
$$

we have

$$
\phi'(0)=
\nabla f(x_k)^Td_k
$$

and

$$
\phi''(0)=
d_k^TH(x_k)d_k.
$$

Therefore,

$$
\boxed{
\phi(\alpha)
\approx
f(x_k)
+
\alpha\nabla f(x_k)^Td_k
+
\frac12\alpha^2d_k^TH(x_k)d_k
}
$$

This equation is extremely important.

It tells us what happens when we move along the chosen direction.

---

# 10. Now choose the gradient descent direction

For gradient descent,

$$
d_k=-\nabla f(x_k).
$$

Substitute this into the Taylor approximation.

The first-order term becomes

$$
\nabla f(x_k)^Td_k=
-\|\nabla f(x_k)\|^2.
$$

Therefore,

$$
\boxed{
\phi(\alpha)
\approx
f(x_k)
-
\alpha\|\nabla f(x_k)\|^2
+
\frac12\alpha^2
\nabla f(x_k)^T
H(x_k)
\nabla f(x_k)
}
$$

Now we can see exactly what the step length does.

---

# 11. Understanding the step length

The first-order term is

$$
-\alpha\|\nabla f(x_k)\|^2.
$$

It is negative.

So initially, increasing $\alpha$ decreases the function.

But the second-order term is

$$
\frac12\alpha^2
\nabla f(x_k)^TH(x_k)\nabla f(x_k).
$$

For positive curvature, this term is positive.

So there is a competition:

$$
\boxed{
\text{first-order decrease}
\quad\leftrightarrow\quad
\text{second-order curvature}
}
$$

For very small $\alpha$, the first-order term dominates.

For larger $\alpha$, curvature becomes increasingly important.

This explains why **the step length matters**.

---

# 12. The line-search picture

Imagine fixing the current point $x_k$ and direction $d_k$.

We obtain a one-dimensional curve:

$$
\phi(\alpha)=f(x_k+\alpha d_k).
$$

Graphically:

```text
f
│
│       ●
│      / \
│     /   \
│    /     \
│   /       \
│  ●
│
└──────────────── α
  0       α*
```

At

$$
\alpha=0
$$

we are at the current point $x_k$.

Moving initially in the chosen descent direction makes the function decrease.

The best step along that line is

$$
\boxed{
\alpha^* =
\arg\min_{\alpha\geq0}
f(x_k+\alpha d_k).
}
$$

Then

$$
\boxed{
x_{k+1}=x_k+\alpha^*d_k.
}
$$

This is **exact line search**.

---

# 13. Gradient descent and line search are not two unrelated ideas

This distinction is important.

Gradient descent answers:

> **Which direction should I move?**

$$
\boxed{
d_k=-\nabla f(x_k)
}
$$

Line search answers:

> **How far should I move in that direction?**

$$
\boxed{
\alpha_k =
\arg\min_{\alpha\geq0}
f(x_k+\alpha d_k)
}
$$

Together:

$$
\boxed{
x_{k+1} =
x_k+\alpha_kd_k
}
$$

with

$$
\boxed{
d_k=-\nabla f(x_k).
}
$$

Therefore,

$$
\boxed{
x_{k+1} =
x_k-\alpha_k\nabla f(x_k).
}
$$

This is gradient descent with a line-search-selected step size.

---

# 14. Fixed learning rate versus line search

In machine learning, we often write

$$
\boxed{
x_{k+1} =
x_k-\alpha\nabla f(x_k)
}
$$

with a fixed or scheduled $\alpha$.

For example:

$$
\alpha=0.01.
$$

The same step size is used repeatedly.

This is different from exact line search:

$$
\boxed{
\alpha_k =
\arg\min_{\alpha\geq0}
f(x_k-\alpha\nabla f(x_k)).
}
$$

Here $\alpha_k$ is determined from the function at the current iteration.

So:

| Method                                    | Direction           | Step length                         |
| ----------------------------------------- | ------------------- | ----------------------------------- |
| Gradient descent with fixed learning rate | $-\nabla f(x_k)$  | predetermined $\alpha$            |
| Gradient descent + line search            | $-\nabla f(x_k)$  | determined by line search           |
| Newton's method                           | $-H^{-1}\nabla f$ | usually incorporated into direction |

---

# 15. A deeper geometric picture

We can now connect everything we have studied.

Start with the point

$$
x_k.
$$

Choose a direction

$$
d_k.
$$

This creates the line

$$
\boxed{
x_k+\alpha d_k.
}
$$

Restrict the original function to this line:

$$
\boxed{
\phi(\alpha)=f(x_k+\alpha d_k).
}
$$

Then:

$$
\phi'(0) =
\nabla f(x_k)^Td_k
$$

is the **slope along the line**.

And

$$
\phi''(0) =
d_k^TH(x_k)d_k
$$

is the **curvature along the line**.

Therefore:

$$
\boxed{
\text{gradient}
\rightarrow
\text{slope in a direction}
}
$$

and

$$
\boxed{
\text{Hessian}
\rightarrow
\text{curvature in a direction}.
}
$$

Gradient descent chooses the direction

$$
d_k=-\nabla f(x_k)
$$

because this gives the steepest local decrease.

Line search then examines

$$
f(x_k+\alpha d_k)
$$

to determine how far to move.

---

# 16. The complete conceptual chain

The entire idea can be presented to students as one continuous story:

$$
\boxed{
\text{Taylor approximation}
}
$$

gives

$$
f(x+s)
\approx
f(x)+\nabla f(x)^Ts
+\frac12s^TH(x)s.
$$

Then ask:

### Question 1: Which direction decreases $f$?

Look at the first-order term:

$$
\nabla f(x)^Ts.
$$

The best direction is

$$
\boxed{
d=-\nabla f(x).
}
$$

This gives **gradient descent direction**.

---

### Question 2: How far should we move?

Once the direction is fixed, define

$$
\phi(\alpha)=f(x+\alpha d).
$$

Now solve the one-dimensional problem

$$
\boxed{
\min_{\alpha\geq0}\phi(\alpha).
}
$$

This gives **line search**.

---

### Question 3: What controls the behavior along the line?

Taylor gives

$$
\boxed{
\phi(\alpha)
\approx
\phi(0)
+
\alpha\phi'(0)
+
\frac12\alpha^2\phi''(0)
}
$$

where

$$
\boxed{
\phi'(0)=\nabla f(x)^Td
}
$$

and

$$
\boxed{
\phi''(0)=d^TH(x)d.
}
$$

So the gradient determines the **initial slope**, while the Hessian determines the **curvature**.

---

# 17. The key message for students

The important thing is not to memorize

$$
x_{k+1}=x_k-\alpha\nabla f(x_k).
$$

Instead, understand where it comes from:

$$
\boxed{
\text{Taylor approximation}
\rightarrow
\text{local slope}
\rightarrow
\text{descent direction}
\rightarrow
\text{gradient descent}
}
$$

and

$$
\boxed{
\text{chosen direction}
\rightarrow
\text{line}
\rightarrow
\text{1D function}
\rightarrow
\text{line search}
\rightarrow
\text{step length}.
}
$$

The central picture is therefore:

$$
\boxed{
x_k
\;\xrightarrow{\quad d_k=-\nabla f(x_k)\quad}\;
x_k+\alpha d_k
}
$$

and along that line,

$$
\boxed{
\phi(\alpha)=f(x_k+\alpha d_k).
}
$$

Gradient descent chooses **the direction**.

Line search chooses **the distance**.

Taylor approximation explains **why both choices make sense locally**.
