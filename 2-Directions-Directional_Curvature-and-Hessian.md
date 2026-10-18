# Directions, Directional Curvature and the Hessian: A View for Optimization

**Setting.** $f:\mathbb{R}^n\to\mathbb{R}$, $n\ge1$, is twice continuously differentiable. Then the Hessian $H=\nabla^2 f$ is symmetric.

<!-- **Aim.** Explain why one chooses a direction, why the displacement then becomes linear, and how this turns the gradient into "slope in every direction" and the Hessian into "curvature in every direction". The conclusion is convexity read directly from the Hessian. -->

---

## 1. Unlimited directions, one choice

Standing at a point $p\in\mathbb{R}^n$, you can move in infinitely many directions. For $n=2$ they form a circle of directions, for $n=3$ a sphere, and in general the unit sphere $S^{n-1}$. The function $f$ responds differently to each of them.

A derivative measures a **rate of change**, and a rate needs one parameter to change along. With $n$ coordinates free to vary, "the rate of change of $f$" has no single meaning. So, to ask a derivative question, we **choose one direction $v$**.

Once $v$ is chosen, the moving point is no longer free in $\mathbb{R}^n$. It is described by one number $t$, the amount of displacement along $v$ defines a straight line:

$$
\boxed{x(t)=p+t\,v,\qquad t\in\mathbb{R}.}
$$

Changing $v$ rotates the line around $p$, and replacing $v$ by $cv$ keeps the same line but traverses it at a different speed.

**In iterative optimization**, most algorithms update the current point according to the general rule $x_{k+1}=x_k+\alpha_k d_k$, where $d_k$ is the chosen search direction and $\alpha_k$ is the step size.

<!-- ### Why a straight line, and not some curve?

Four reasons, none of which requires an algorithm.

1. **A single finite move is always linear.** Take any other point $q$. Then $q=p+(q-p)$, so $q$ lies on the line with direction $v=q-p$ at $t=1$. Choosing a line never prevents us from reaching a point.
2. **A line is the simplest one-parameter family.** Along it, the displacement is linear in $t$, so every change in $f$ comes from $f$ itself and not from the shape of the path.
3. **Optimization moves linearly.** A step from the current point $x_k$ along a direction $v_k$ with step length $t$ is $x_k+t\,v_k$. The choice of direction is made by the person or the method, not by $f$. Each step defines its own line, through a new base point.
4. **Convexity is defined on lines.** A function is convex if its restriction to every line is convex (Section 8). This is the language the subject uses.

**Two clarifications.**

- The line is not needed to *evaluate* $f$. The value $f(q)$ is defined for every point. The direction is needed to *differentiate*, meaning to reduce an $n$-dimensional question to a scalar one.
- Moving along straight segments is the standard model for unconstrained steps. It is a modelling choice and not a law: gradient flow, constrained problems and manifold methods follow curved paths. Within a small enough neighborhood, any smooth path is well approximated by its tangent line, which is why the straight model is the first one to understand.

--- -->

## 2. Restrict $f$ to the line

Along the chosen line, $f$ becomes a function of one variable:

$$
\boxed{g_v(t)=f(p+t\,v).}
$$

$$
\text{function of } n \text{ variables}\ \longrightarrow\ \text{function of one variable along the chosen direction.}
$$

Now all one-variable tools apply: slope, rate of change of slope, tangent line, parabola.

<!-- ### The family of restrictions is $f$

Different $v$ give different one-variable functions $g_v$. Together they contain exactly the information of $f$, since

$$
f(p+d)=g_d(1)\qquad\text{for every displacement } d .
$$

So the infinite set of choices of $v$ is equivalent to knowing $f$ near $p$, seen one line at a time. The remarkable fact, derived next, is that the **slopes** and **second derivatives** of all these infinitely many functions are determined by a finite number of values at $p$.

Note also the scaling rule $g_{cv}(t)=g_v(ct)$. Doubling $v$ does not change the line, it only runs along it twice as fast. Therefore $g_{cv}'(0)=c\,g_v'(0)$ and $g_{cv}''(0)=c^2g_v''(0)$. This already hints that the slope is linear and the second derivative quadratic in $v$. -->

<!-- --- -->

## 3. First derivative: slope along $v$ is linear in $v$

Assume $f$ is differentiable at $p$. By the chain rule,

Starting from
$$
g_v(t)=f(p+tv)
$$

we want to calculate how $f$ changes as we move along the line.

**Step 1: Express each coordinate**

Since $x,p,v\in\mathbb{R}^n$, write

$$
x=(x_1, \ldots, x_n),\qquad p=(p_1,\ldots,p_n),\qquad v=(v_1,\ldots,v_n)
$$

Then, writing $x$ as function of $t$ gives:

$$
x(t)=p+tv=(p_1+tv_1,\ldots,p_n+tv_n).
$$

Therefore,

$$
g_v(t)=f(p_1+tv_1,\ldots,p_n+tv_n).
$$

Key idea: $f$ depends on $n$ coordinates, but along this line, every coordinate depends on the single variable $t$.

Step 2: Apply the chain rule

Differentiate $g_v(t)$ with respect to $t$. Since each coordinate changes with $t$, the chain rule gives

$$
\frac{dg_v}{dt} =
\underbrace{\frac{\partial f}{\partial x_1}\frac{dx_1}{dt}}_{\text{contribution from }x_1} + \ldots +
\underbrace{\frac{\partial f}{\partial x_n}\frac{dx_n}{dt}}_{\text{contribution from }x_n}.
$$

Since $x(t)=p+tv$ and $x_i(t)=p_i+tv_i$, we obtain

$$
g_v'(t)=\sum_{i=1}^{n}
\frac{\partial f}{\partial x_i}(p+tv)
\frac{d}{dt}(p_i+tv_i).
$$

Now,

$$
\frac{d}{dt}(p_i+tv_i)=v_i,
$$

because $p_i$ and $v_i$ are constants.
Substituting yields

$$
g_v'(t)=\sum_{i=1}^n \frac{\partial f}{\partial x_i}(p+tv)\,v_i
$$

Each partial derivative measures the sensitivity of $f$ to one coordinate, while $v_i$ measures how quickly that coordinate changes along the line. The total rate of change is the sum of these contributions.

Step 3: Express the sum using the gradient

Recall that the gradient is the vector of partial derivatives:

$$
\nabla f(x)=
\begin{pmatrix}
\frac{\partial f}{\partial x_1}(x)\\\\
\vdots\\\\
\frac{\partial f}{\partial x_n}(x)
\end{pmatrix},
\qquad
v=\begin{pmatrix}v_1\\\\\vdots\\\\v_n\end{pmatrix}.
$$

Their scalar product is exactly the sum obtained above. Hence,

$$
\boxed{g_v'(t)=\nabla f(p+tv)^T v.}
$$

Step 4: Evaluate at the starting point:

At $t=0$, we are at the initial point:

$$
p+0v=p.
$$

Therefore,

$$
\boxed{g_v'(0)=\nabla f(p)^T v.}
$$

Interpretation: $g_v'(0)$ calculate the directional derivative at the posiiton $p$:which is the gradient projected onto the chosen direction. It tells us the instantaneous rate at which $f$ changes when we leave $p$ along the direction $v$.
If $v$ is a unit vector, $g_v'(0)$ is the rate of change per unit of distance traveled. If $v$ is not a unit vector, the rate is also scaled by its magnitude.


<!-- Reading:

- $\nabla f(p)^Tv$ is the **rate of change of $f$ at $p$ per unit of $t$** along $v$. For $\|v\|=1$ it is the directional derivative.
- It is **linear in $v$**. Knowing the $n$ numbers $\partial f/\partial x_i(p)$ gives the slope in **every** direction at once. This is the content of the gradient: $n$ numbers summarizing infinitely many slopes.
- Among unit directions the slope is largest for $v=\nabla f/\|\nabla f\|$, with value $\|\nabla f\|$, and the slope is zero for directions perpendicular to $\nabla f$.

*Remark on hypotheses.* Having a derivative along every line does not by itself guarantee linearity in $v$. The formula $g_v'(0)=\nabla f(p)^Tv$ needs $f$ to be differentiable at $p$. With continuous partial derivatives this is automatic. -->

---

## 4. Second derivative: drift of the slope along $v$ is quadratic in $v$

We want to differentiate the first derivative once more to understand how the rate of change itself changes as we move along the line.

$$
g_v'(t)=\sum_i f_i(p+tv)\,v_i
$$ 
<!-- once more, where $f_i=\partial f/\partial x_i$. By the chain rule each $f_i$ changes along the line at rate $\sum_j f_{ij}\,v_j$, so -->

**Step 1: Identify what changes with $t$**

In each term,

$$
f_i(p+tv)\,v_i,
$$

the factor $v_i$ is constant, because the direction $v$ is fixed. However, $f_i(p+tv)$ changes because we move to different points along the line.
Therefore, differentiating gives

$$
g_v''(t)=\sum_{i=1}^{n}
\frac{d}{dt}f_i(p+tv)\,v_i.
$$

The question now becomes: how do we differentiate $f_i(p+tv)$?

**Step 2: Apply the chain rule again**

Recall that $f_i$ is a function of $n$ coordinates:

$$
f_i=f_i(x_1,\ldots,x_n).
$$

As we move along the line, every coordinate changes according to $x_j(t)=p_j+tv_j$.

The chain rule therefore gives

$$
\frac{d}{dt}f_i(p+tv)=
\sum_{j=1}^{n}
\frac{\partial f_i}{\partial x_j}(p+tv)
\frac{dx_j}{dt}.
$$

Since

$$
\frac{\partial f_i}{\partial x_j}=
\frac{\partial^2 f}{\partial x_j\partial x_i}
=f_{ij},
\qquad
\frac{dx_j}{dt}=v_j,
$$

we obtain

$$
\boxed{
\frac{d}{dt}f_i(p+tv) =
\sum_{j=1}^{n} f_{ij}(p+tv)v_j
}
$$

Interpretation: The rate of change of the $i$-th partial derivative is the sum of its sensitivities to every coordinate, weighted by how quickly those coordinates move.

**Step 3: Substitute into the second derivative**

Returning to

$$
g_v''(t)=\sum_{i=1}^{n}
\frac{d}{dt}f_i(p+tv)\,v_i,
$$

Substitute the expression from Step 2:

$$
g_v''(t) =
\sum_{i=1}^{n}
\left(\sum_{j=1}^{n}f_{ij}(p+tv)v_j\right)v_i.
$$

Rearranging the factors 

$$
\boxed{
g_v''(t)=
\sum_{i=1}^{n}\sum_{j=1}^{n}
v_i f_{ij}(p+tv)v_j
}
$$

**Step 4: Recognize the structure**

The second derivative involves three elements:

$$
\underbrace{v_i}_{\text{first direction component}}
\quad
\underbrace{f_{ij}}_{\text{second derivative}}
\quad
\underbrace{v_j}_{\text{second direction component}}.
$$

The second derivative $f_{ij}$ describes how the slope in coordinate $i$ changes when we move in coordinate $j$. The two direction components account for how our chosen movement interacts with that change.

Collecting all second derivatives into the Hessian matrix $H$, we can write the double sum more compactly as

$$
\boxed{
g_v''(t)=v^T H(p+tv)v
}
$$

At the starting point $t=0$,

$$
\boxed{g_v''(0)=v^TH(p)v.}
$$

The central idea: The first derivative measures the slope along the chosen direction. Differentiating it again measures how that slope changes along the same direction. The Hessian captures all the changes in coordinate slopes, and the quadratic form $v^THv$ combines them into the curvature along the chosen line.

**Vocabulary**

- $v^THv$ is the **second directional derivative** of $f$ along $v$. It scales with $\|v\|^2$.
- For a **unit** vector $u$, $u^THu$ is the second derivative per unit length, and only then is it called the **directional curvature**.




<!-- $$
g_v''(t)=\sum_{i=1}^n\sum_{j=1}^n v_i\,f_{ij}(p+tv)\,v_j .
$$

Collecting the second partial derivatives in the Hessian matrix $H_{ij}=f_{ij}$:

$$
\boxed{g_v''(t)=v^T H(p+tv)\,v,\qquad g_v''(0)=v^T H(p)\,v.}
$$

Reading:

- The vector $Hv$ is the **rate of change of the gradient** per unit of $t$ along the line: $\frac{d}{dt}\nabla f(p+tv)\big|_{t=0}=H(p)\,v$.
- Multiplying by $v^T$ gives the **drift of the slope along $v$**. The part of $Hv$ perpendicular to $v$ is the gradient turning sideways and does not enter.
- $v^THv$ is **quadratic in $v$**: a symmetric $n\times n$ matrix has $\tfrac{n(n+1)}2$ independent entries, and these determine the second derivative in every direction.

| Quantity | Number of values needed | Dependence on $v$ |
|---|---|---|
| slope $g_v'(0)$ | $n$ (gradient) | linear |
| second derivative $g_v''(0)$ | $\tfrac{n(n+1)}2$ (Hessian) | quadratic | -->

**Example**

For $n=2$ the expansion reads

$$
v^THv=f_{xx}v_1^2+2f_{xy}v_1v_2+f_{yy}v_2^2 ,
$$

with the cross term appearing twice because $H_{12}=H_{21}$.


---

## 5. What happens inside $v^THv$

The key idea is that $v^THv$ is calculated in two steps: first, determine how the gradient changes when we move along $v$; then, determine how much of that change affects the slope in the same direction $v$.

**1. Start with the gradient**

$$
v\ \xrightarrow{\ H\ }\ Hv\ \xrightarrow{\ v^T\ }\ v^THv .
$$

At a point $p$, the gradient is

$$
\nabla f(p)=
\begin{pmatrix}
f_{x_1}\\\\
\vdots\\\\
f_{x_n}
\end{pmatrix}.
$$

It contains the slopes of $f$ along all coordinate directions.
When we move a small distance in direction $v$, the gradient changes. The Hessian describes this change:

$$
d(\nabla f)=H\,d\mathbf X.
$$

If the displacement is $d\mathbf X=v\,dt$, then

$$
\frac{d}{dt}\nabla f(p+tv)\bigg|_{t=0}
=Hv.
$$

Thus, $Hv$ is the rate at which the entire gradient vector changes as we move along $v$.

**2. Why multiply $Hv$ by $v^T$?**

The vector $Hv$ tells us how every component of the gradient changes. But our question is more specific:
How does the slope of the function along our chosen direction $v$ change?
The slope along $v$ is

$$
g_v'(t)=\nabla f(p+tv)^T v.
$$

Differentiate once more. Since $v$ is fixed,

$$
g_v''(0) =
\left(\frac{d}{dt}\nabla f(p+tv)\bigg|_{t=0}\right)^T v.
$$

Substitute $Hv$:

$$
g_v''(0)=(Hv)^Tv.
$$

For a symmetric Hessian, this equals

$$
\boxed{g_v''(0)=v^THv.}
$$

The multiplication by $v^T$ takes the gradient's change and measures its component along the direction $v$. That is exactly the information needed to determine how the slope along $v$ changes.

<!-- | Quantity | Type | Meaning |
|---|---|---|
| $Hv$ | vector | how the gradient changes per unit step along $v$ |
| $v^T(Hv)$ | scalar | how the slope **along $v$** changes, which is $g_v''(0)$ |

The first operation is the Hessian acting on the direction. The second keeps only the component along the direction chosen, which is what a one-variable second derivative can see. -->

---

<!-- ## 6. Examples

### Example 1: diagonal Hessian, $f(x,y)=x^2+2y^2$

$$
\nabla f=\begin{pmatrix}2x\\4y\end{pmatrix},\qquad
H=\begin{pmatrix}2&0\\0&4\end{pmatrix},\qquad
v^THv=2v_1^2+4v_2^2 .
$$

| Direction $v$ | $v^THv$ | Direct check along the line from the origin |
|---|---|---|
| $(1,0)$ | $2$ | $f(t,0)=t^2$, so $g''=2$ |
| $(0,1)$ | $4$ | $f(0,t)=2t^2$, so $g''=4$ |
| $(1,1)$ | $6$ | $f(t,t)=3t^2$, so $g''=6$ |
| $(1,1)/\sqrt2$ | $3$ | $f(t/\sqrt2,t/\sqrt2)=\tfrac32t^2$, so $g''=3$ |

The last two rows show the scaling with $\|v\|^2$: the same line with a unit vector gives half of $6$.

### Example 2: non-diagonal Hessian, $f(x,y)=x^2+3xy$

$$
\nabla f=\begin{pmatrix}2x+3y\\3x\end{pmatrix},\qquad
H=\begin{pmatrix}2&3\\3&0\end{pmatrix},\qquad
v^THv=2v_1^2+6v_1v_2 .
$$

| Direction $v$ | $v^THv$ | Direct check |
|---|---|---|
| $(1,0)$ | $2$ | $f(t,0)=t^2$ |
| $(0,1)$ | $0$ | $f(0,t)=0$ |
| $(1,1)$ | $8$ | $f(t,t)=4t^2$ |
| $(1,-1)$ | $-4$ | $f(t,-t)=-2t^2$ |

The two coordinate axes show curvature $2$ and $0$, which suggests "no downward bending". The diagonal $(1,-1)$ bends **downward**. The off-diagonal entry is exactly what the axes cannot see.

The extreme values of $v^THv$ over unit vectors are the eigenvalues of $H$ (Rayleigh bound $\lambda_{\min}\|v\|^2\le v^THv\le\lambda_{\max}\|v\|^2$):

$$
\lambda=1\pm\sqrt{10}\approx 4.16,\ -2.16 .
$$

The negative eigenvalue confirms a downward direction. Since $\nabla f(0,0)=0$, the origin is a **saddle**.

---

## 7. Definiteness as curvature in every direction

For symmetric $H$, using $g_v''(0)=v^THv$:

| Condition | Meaning for $f$ at $p$ | Eigenvalues |
|---|---|---|
| $v^THv\ge0$ for all $v$ ($H\succeq0$) | no direction bends downward | all $\ge0$ |
| $v^THv>0$ for all $v\ne0$ ($H\succ0$) | strictly upward in every direction, no flat direction | all $>0$ |
| $v^THv<0$ for some $v$ | some direction bends downward | some $<0$ |
| both signs occur ($H$ indefinite) | up in some directions, down in others | mixed signs |

**Saddle.** For $f=x^2-y^2$, $H=\mathrm{diag}(2,-2)$ and $v^THv=2v_1^2-2v_2^2$. Along $x$ the restriction is $t^2$ (up), along $y$ it is $-t^2$ (down). At a **critical point** an indefinite Hessian means a saddle. Away from a critical point the slope term also matters, so indefiniteness alone does not decide.

**Flat direction.** For $f=e^{x+2y}$, $H=e^{x+2y}\begin{pmatrix}1&2\\2&4\end{pmatrix}$, which is PSD of rank one. For $v=(2,-1)$ the restriction is constant, and $v^THv=0$: zero curvature along that direction, yet no direction bends downward.

--- -->

