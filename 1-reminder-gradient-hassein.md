# Reminder: derivative, second derivative, gradient, Hessian

## 1. One variable

**First derivative = slope.**

$$
f'(x)=\lim_{\Delta x\to0}\frac{\Delta f}{\Delta x}
$$

It is the rate at which $f$ changes per unit of $x$. Its sign tells you whether $f$ rises or falls.

**First differential.** For a chosen displacement $dx$,

$$
df=f'(x)\,dx .
$$

Read it as "slope × displacement". $df$ is exactly linear in $dx$ and has no second-order part. The notation $\frac{df}{dx}$ is the coefficient linking the two.

**Second derivative = rate of change of the slope.**

$$
f''(x)=\frac{d}{dx}\big(f'(x)\big)
$$

- $f''>0$: the slope increases as you move right (the curve bends upward).
- $f''<0$: the slope decreases (the curve bends downward).
- $f''=0$: the slope is momentarily not changing.

**Second differential.** Differentiate $df=f'(x)\,dx$ again, holding $dx$ fixed:

$$
d(f')=f''\,dx,\qquad d^2f=dx\,\big(f''dx\big)=f''(x)\,dx^2 .
$$

It is the change in slope, $f''dx$, multiplied once more by the displacement $dx$. It has the same shape as the first order, "slope × displacement", but with the slope replaced by its change.

**Combining the two orders.** The slope is not constant during the move. It starts at $f'$ and drifts linearly at rate $f''$:

$$
\text{slope at distance } s \;=\; f' + f''\,s .
$$

The real change $\Delta f$ is the area under this slope line, a rectangle plus a triangle:

$$
\Delta f=\underbrace{f'\,dx}_{df:\ \text{rectangle}}+\underbrace{\tfrac12 f''\,dx^2}_{\tfrac12 d^2f:\ \text{triangle}}
$$

- $df$ assumes the initial slope holds for the whole move.
- $\tfrac12 d^2f$ corrects for the slope drift.
- The $\tfrac12$ is the triangle's area. The drift grows from $0$ to $f''dx$, so it is not at full size for the whole move.

Equivalently, $\Delta f$ is the average slope times $dx$, and for a linearly drifting slope the average is the midpoint value $f'+\tfrac12 f''dx$.

**Check.** For $f=x^2$ at $x=2$ with $dx=0.1$: $df=0.4$ and $\tfrac12 d^2f=\tfrac12(2)(0.01)=0.01$, so $\Delta f=0.41$, which is exact because $f''$ is constant.

---

## 2. Two variables

Now $f(x,y)$ has one slope per direction:

$$
f_x=\frac{\partial f}{\partial x},\qquad f_y=\frac{\partial f}{\partial y}.
$$

**First differential.** Each direction contributes slope × displacement, and the contributions add:

$$
df=f_x\,dx+f_y\,dy .
$$

**Second differential, step by step.** Differentiate $df$ with $dx,dy$ held fixed. Only the slopes vary, and both can change with position:

$$
d^2f=df_x\,dx+df_y\,dy,
\qquad
df_x=f_{xx}\,dx+f_{xy}\,dy,\quad df_y=f_{yx}\,dx+f_{yy}\,dy .
$$

Substituting gives four terms:

$$
d^2f=f_{xx}\,dx\,dx+f_{xy}\,dy\,dx+f_{yx}\,dx\,dy+f_{yy}\,dy\,dy .
$$

Each term has the pattern (second derivative) × (displacement in one variable) × (displacement in another). With continuous second derivatives, $f_{xy}=f_{yx}$, so the mixed terms merge:

$$
d^2f=f_{xx}\,dx^2+2f_{xy}\,dx\,dy+f_{yy}\,dy^2 .
$$

Reading the terms:

- $f_{xx}$: the $x$-slope changes when you move in $x$.
- $f_{yy}$: the $y$-slope changes when you move in $y$.
- $f_{xy}$: the $x$-slope changes when you move in $y$, and the $y$-slope changes when you move in $x$. This is the interaction term, and it counts twice.

**Combining the orders.** Walk the straight path $(x+t\,dx,\;y+t\,dy)$ for $t\in[0,1]$ and set $\varphi(t)=f(x+t\,dx,\,y+t\,dy)$. This is a one-variable function, so the one-variable picture applies unchanged:

- $\varphi'(0)=df$ is the initial slope along the path.
- $\varphi''(0)=d^2f$ is the drift of that slope.
- The rectangle plus triangle argument gives

$$
\Delta f\approx df+\tfrac12\,d^2f .
$$

---

## 3. Gradient and Hessian: the same objects in vector form

Let $d\mathbf X=(dx,dy)^T$.

**Gradient: the coefficients of $df$.**

$$
df=\begin{pmatrix}f_x&f_y\end{pmatrix}\begin{pmatrix}dx\\\\dy\end{pmatrix}=\nabla f^{\,T}d\mathbf X,
\qquad
\nabla f=\begin{pmatrix}f_x\\\\f_y\end{pmatrix}.
$$

It turns a displacement into a change in $f$. It points in the direction of steepest increase, and its norm is the steepest slope.

**Hessian: the coefficients of $d(\nabla f)$.** The gradient is the vector of slopes, so differentiating it means differentiating each slope. Each slope depends on both coordinates, so each one changes with both $dx$ and $dy$.

*Step 1: differentiate each component.*

$$
d(\nabla f)=\begin{pmatrix}df_x\\\\ df_y\end{pmatrix}
=\begin{pmatrix}f_{xx}\,dx+f_{xy}\,dy\\\\[2pt] f_{yx}\,dx+f_{yy}\,dy\end{pmatrix}.
$$

*Step 2: recognise this column as a matrix times a vector.* Each row is a row of second derivatives dotted with $(dx,dy)^T$:

$$
\begin{pmatrix}f_{xx}\,dx+f_{xy}\,dy\\\\[2pt] f_{yx}\,dx+f_{yy}\,dy\end{pmatrix}
=\begin{pmatrix}f_{xx}&f_{xy}\\\\ f_{yx}&f_{yy}\end{pmatrix}\begin{pmatrix}dx\\\\dy\end{pmatrix}.
$$

*Step 3: name the matrix.* It is the Hessian:

$$
\boxed{d(\nabla f)=H\,d\mathbf X},\qquad H=\begin{pmatrix}f_{xx}&f_{xy}\\\\ f_{yx}&f_{yy}\end{pmatrix}.
$$

The Hessian turns a displacement into a change in the gradient. This is a **vector** with one entry per slope. No row $d\mathbf X^T$ appears yet, because nothing has been multiplied by a displacement.

**Second differential as a quadratic form.** The row $d\mathbf X^T$ enters when you differentiate the scalar $df=\nabla f^Td\mathbf X$. With $d\mathbf X$ held fixed, only $\nabla f$ varies:

$$
d^2f=d(\nabla f)^T\,d\mathbf X=(H\,d\mathbf X)^Td\mathbf X=d\mathbf X^T H^T d\mathbf X=d\mathbf X^T H\,d\mathbf X,
$$

using $H=H^T$. Written out in the form of the earlier sum:

$$
d^2f=dx\,(f_{xx}dx+f_{xy}dy)+dy\,(f_{yx}dx+f_{yy}dy)
$$

$$
\boxed{d^2f=d\mathbf X^T H\,d\mathbf X}
$$

Expanding returns the four-term sum, then $f_{xx}dx^2+2f_{xy}dxdy+f_{yy}dy^2$.

**Do not mix the two objects:**

| Object | Type | Its differential | Type |
|---|---|---|---|
| $\nabla f$ | vector | $d(\nabla f)=H\,d\mathbf X$ | vector: how the slopes change |
| $df=\nabla f^Td\mathbf X$ | scalar | $d(df)=d^2f=d\mathbf X^TH\,d\mathbf X$ | scalar: how $df$ changes |

The scalar $d^2f$ is the change in the gradient projected onto the displacement.

**The structure mirrors the first order:**

- $df=\nabla f^T d\mathbf X$: gradient row times displacement, which is slope times displacement.
- $d^2f=d\mathbf X^T(H\,d\mathbf X)$: displacement row times the change in the gradient, which is the same pattern with the slope replaced by its drift.
- In one variable, $d^2f=dx\,(f''dx)=f''dx^2$, and $H$ collapses to the scalar $f''$.

**What the matrix is bookkeeping.** Entry $H_{ij}=f_{x_ix_j}$ is the rate at which slope $i$ changes when you move in direction $j$. The quadratic form sums, over every pair $(i,j)$,

$$
\underbrace{dx_i}_{\text{displacement}}\times\underbrace{H_{ij}}_{\text{slope } i \text{ drift per unit } x_j}\times\underbrace{dx_j}_{\text{displacement}} .
$$

The off-diagonal term appears twice, once as $(i,j)=(1,2)$ and once as $(2,1)$, which is where the factor 2 in $2f_{xy}\,dx\,dy$ comes from.

**Dictionary:**

| One variable | Several variables |
|---|---|
| $f'$ | $\nabla f$ |
| $f''$ | $H$ |
| $df=f'\,dx$ | $df=\nabla f^T d\mathbf X$ |
| $d(f')=f''\,dx$ | $d(\nabla f)=H\,d\mathbf X$ |
| $d^2f=f''dx^2$ | $d^2f=d\mathbf X^T H\,d\mathbf X$ |
| $\Delta f\approx f'dx+\tfrac12f''dx^2$ | $\Delta f\approx\nabla f^Td\mathbf X+\tfrac12d\mathbf X^TH\,d\mathbf X$ |

The scalar product $f''dx^2$ becomes the quadratic form $d\mathbf X^THd\mathbf X$.

---

## 4. Meaning to remember

- **$\nabla f^Td\mathbf X$** is the predicted change if the slopes stayed frozen.
- **$\tfrac12 d\mathbf X^TH\,d\mathbf X$** is the correction due to the slopes changing during the move.
- For a unit direction $u$, $u^THu$ is the second derivative of $f$ along that line, meaning how sharply $f$ bends in that direction.
- $d\mathbf X^THd\mathbf X>0$ for all $d\mathbf X\neq0$ means $f$ bends upward in every direction, which is a bowl.
- $H$ is symmetric when the second derivatives are continuous.

---

## 5. Numerical check

$f(x,y)=x^2+3xy$ at $(1,1)$ with $d\mathbf X=(0.1,0.2)^T$.

- $\nabla f=(2x+3y,\;3x)=(5,3)$, so $df=0.5+0.6=1.1$.
- $H=\begin{pmatrix}2&3\\\\3&0\end{pmatrix}$, so

$$
d(\nabla f)=H\,d\mathbf X=\begin{pmatrix}0.2+0.6\\\\0.3+0\end{pmatrix}=\begin{pmatrix}0.8\\\\0.3\end{pmatrix}\ (\text{vector}),
\qquad
d^2f=d\mathbf X^T(H\,d\mathbf X)=0.1(0.8)+0.2(0.3)=0.14\ (\text{scalar}).
$$

- Direct check of the vector: $df_x=2(0.1)+3(0.2)=0.8$ and $df_y=3(0.1)=0.3$.
- The expanded form agrees: $2(0.01)+2\cdot3(0.02)+0=0.14$.
- $\Delta f=1.1+\tfrac12(0.14)=1.17$, which matches the exact $f(1.1,1.2)-f(1,1)=1.17$. It is exact here because $f$ is quadratic, so $H$ is constant.

---

<!-- ## 6. Caveats

1. **Independent displacement.** This treats $dx$ as an independent displacement. If $x$ itself depends on another variable, $d^2f$ gains an extra term $f'\,d^2x$.
2. **Where the formula holds.** $\Delta f=df+\tfrac12 d^2f$ is exact only when the second derivatives are constant along the path. Otherwise it is a local approximation, and any statement about how large the error is belongs to Taylor territory.
3. **Relation to Taylor.** The formula for $\Delta f$ is the second-order Taylor formula. The route here is slope, drift of slope, and area, so no series is needed to derive or interpret it. -->
