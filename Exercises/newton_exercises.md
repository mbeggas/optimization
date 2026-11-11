# Newton's Method — Exercises 

**Scope:** Part I (Newton's method for root-finding) and Part II (Newton's method for optimization).

**Objective:** Understand how Newton's linearization becomes a method for optimization, how the Hessian modifies the gradient direction, and why Newton behaves differently from gradient descent.

---

## Exercise 01 — Newton for a System of Nonlinear Equations (4 points)

**Aim:** Understand the Jacobian as the linear transformation predicting changes in a vector-valued function.

Consider the nonlinear system

$$
F(x,y)=\begin{pmatrix}x^2+y^2-4\\\\xy-1\end{pmatrix}=\begin{pmatrix}0\\\\0\end{pmatrix}.
$$

Starting from

$$
\mathbf{x}_0=\begin{pmatrix}2\\\\0.8\end{pmatrix},
$$

answer the following questions.

1. Compute $F(\mathbf{x}_0)$ and the Jacobian $J(\mathbf{x}_0)$.
2. Write the first-order approximation $F(\mathbf{x}_0+p)\approx F(\mathbf{x}_0)+J(\mathbf{x}_0)p$. Explain what the vector $J(\mathbf{x}_0)p$ represents.
3. Solve the Newton system $J(\mathbf{x}_0)p=-F(\mathbf{x}_0)$ and compute $\mathbf{x}_1=\mathbf{x}_0+p$.
4. Compare $\|F(\mathbf{x}_0)\|_2$ and $\|F(\mathbf{x}_1)\|_2$. Has the residual decreased?
5. Explain geometrically why the Newton step is the intersection of the two *linearized* equations $F_1=0$ and $F_2=0$.

**Conceptual question:** If $J(\mathbf{x}_0)$ is singular, why might the Newton correction have no solution or infinitely many solutions?

---


## Exercise 02 — Directional Curvature, Positive Definiteness, and Newton's Method (5 points)

**Aim:** Understand positive definiteness through the curvature of a function in different directions, then explain how the Hessian transforms the gradient and how Newton's method corrects this transformation.

Consider the quadratic objective

$$
f(\mathbf{x})=\frac12\mathbf{x}^TA\mathbf{x},
\qquad
A=\begin{pmatrix}5&4\\\\4&5\end{pmatrix},
$$

with starting point

$$
\mathbf{x}_0=\begin{pmatrix}1\\\\0\end{pmatrix}.
$$

### Part A — Directional curvature and positive definiteness

1. Compute the gradient and Hessian of $f$. Verify that the Hessian is symmetric.

2. Consider an arbitrary nonzero direction $d=(a,b)^T$ and the one-dimensional restriction

   $$
   \phi(t)=f(\mathbf{x}+td).
   $$

   Show that

   $$
   \phi''(0)=d^TAd.
   $$

   Explain the geometric meaning of this quantity.

3. Calculate $d^TAd$ and show that it is strictly positive for every nonzero direction $d$. What does this imply about the curvature of the quadratic surface and the positive definiteness of $A$?

4. Verify the positive definiteness of $A$ directly using the matrix criterion for a symmetric $2\times2$ matrix. Determine the unique minimizer $\mathbf{x}^*$.

### Part B — Principal curvature directions and gradient distortion

5. Find the eigenvalues and an orthonormal eigenvector basis of $A$. Identify the high-curvature and low-curvature directions. Explain how the eigenvalues relate to directional curvature.

6. Express the displacement $e=\mathbf{x}_0-\mathbf{x}^*$ in the eigenvector basis.

7. Compute the gradient $g_0=Ae$ in both the original coordinates and the eigenvector basis. Explain which component is amplified most strongly and why.

### Part C — Gradient descent versus Newton's method

8. Compute the gradient descent direction $d_{\mathrm{GD}}=-g_0$ and the Newton direction $p_N=-A^{-1}g_0$.

9. Calculate the angle between each direction and the direct displacement $\mathbf{x}^*-\mathbf{x}_0$.

10. Show that a full Newton step reaches the minimizer exactly.

11. **Geometric interpretation:** Sketch elliptical level sets of $f$, the two principal curvature directions, $-g_0$, and $p_N$ on the same diagram. Explain why the two directions differ.

**Conceptual question:** Does $A^{-1}$ merely rotate the gradient? Explain its action on the gradient components along the eigenvectors and how this corrects the curvature-dependent distortion.


## Exercise 03 — Why Gradient Descent Zigzags but Newton Does Not (5 points)

**Aim:** Connect Hessian eigenvalues to the convergence behavior of gradient descent.

Consider

$$
f(x,y)=\frac12(20x^2+2y^2),
$$

with starting point $\mathbf{x}_0=(1,1)^T$. Use gradient descent with constant learning rate $\alpha=0.09$.

1. Compute the gradient and Hessian. Write the gradient descent update in matrix form.
2. Find the eigenvalues and eigenvectors of the Hessian. Identify the stiff and soft directions.
3. Derive a recurrence for each coordinate, and write $x_k$ and $y_k$ explicitly as functions of $k$.
4. Calculate the first four iterations $\mathbf{x}_1,\mathbf{x}_2,\mathbf{x}_3,\mathbf{x}_4$.
5. Explain why the high-curvature component alternates in sign while the low-curvature component decreases slowly.
6. Derive the interval of constant learning rates for which gradient descent converges from every starting point.
7. Find the constant learning rate minimizing the worst-case error contraction factor. Compute that factor for this objective.
8. Compute one full Newton step from $\mathbf{x}_0$ and compare it with gradient descent.

**Conceptual question:** Why can a single scalar learning rate generally not eliminate both error components in one iteration, whereas Newton can?

---

## Exercise 04 — Does Newton Always Find a Minimum? (5 points)

**Aim:** Distinguish solving the stationarity equation from minimizing an objective, and interpret singular Hessians.

Consider

$$
\begin{aligned}
f_1(x,y)&=x^2+y^2,\\\\
f_2(x,y)&=-x^2-y^2,\\\\
f_3(x,y)&=x^2-y^2,\\\\
f_4(x,y)&=x^2+y^4.
\end{aligned}
$$

1. Find all stationary points of the four functions.
2. Compute the Hessian at each stationary point.
3. Use Hessian eigenvalues to classify the stationary points whenever the second-order test is conclusive.
4. Starting from $\mathbf{x}_0=(1,1)^T$, compute one Newton step for $f_1$, $f_2$, and $f_3$.
5. Explain why Newton reaches a stationary point in one step for all three functions, although only one has a minimum there.
6. For $f_4$, compute the Hessian at $(0,0)$. What does its zero eigenvalue imply about second-order bending in the corresponding direction?
7. Can you still establish whether $(0,0)$ is a minimum of $f_4$? Justify your answer without relying only on the Hessian test.

**Conceptual question:** Why is positive definiteness of the Hessian important when using Newton's method for optimization rather than for general root-finding?
