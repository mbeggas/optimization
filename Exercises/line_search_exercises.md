# TD: Gradient Descent, Exact Line Search and Armijo Backtracking

**Notation.** For a point $x_k$ and a direction $d_k$, the line-search function is $\phi(\alpha)=f(x_k+\alpha d_k)$, and by the chain rule $\phi'(\alpha)=\nabla f(x_k+\alpha d_k)^T d_k$.

| Exercise | Topic | Function |
|---|---|---|
| 1 | Gradient descent with a fixed step, then exact line search | $f(x,y)=\ln(1+x^2+2y^2)$ |
| 2 | Armijo backtracking (no closed-form exact step) | $f(x,y)=x^4+y^2$ |

---

# Exercise 1: Fixed step, then exact line search

Let
$$
f(x,y)=\ln\!\left(1+x^2+2y^2\right),\qquad p=(1,1).
$$

**Aim.** Show that the descent direction alone is not enough: a fixed step can be good, bad, or even increase $f$, while a step chosen by exact line search adapts at each iteration. Also see the geometric meaning of the condition $\phi'(\alpha^*)=0$.

## Statement

1. Compute $\nabla f(x,y)$, then $\nabla f(p)$ and $d_0=-\nabla f(p)$.
2. **Fixed step $\alpha=1$.** Write the update $x_{k+1}=x_k-\alpha\nabla f(x_k)$. Compute $x_1$ and $x_2$, and $f(x_0),f(x_1),f(x_2)$.
3. **Fixed step $\alpha=2$.** Compute $x_1$ and $x_2$ and the corresponding values of $f$. What do you observe?
4. **Exact line search from $p$.** Write $\phi(\alpha)=f(p+\alpha d_0)$ and simplify it. Check that $\phi'(0)=-\|\nabla f(p)\|^2$. Find $\alpha_0^*$, then $x_1$ and $f(x_1)$.
5. Compute $\nabla f(x_1)^Td_0$. What does it say about the next direction $d_1=-\nabla f(x_1)$?
6. Perform a second exact line search from $x_1$. Give $\alpha_1^*$, $x_2$ and $f(x_2)$.
7. Compare, in a table, the values $f(x_0),f(x_1),f(x_2)$ for the fixed step $\alpha=1$ and for exact line search.
8. *(Remark)* For a pure quadratic, the fixed-step update is linear and the condition on $\alpha$ is a single inequality. Why is this not the case here? What local condition on $\alpha$ appears near the minimizer $(0,0)$?

## Solution

### Question 1: Gradient and direction

Write $q=x^2+2y^2$, so $f=\ln(1+q)$ and
$$
\nabla f(x,y)=\frac{1}{1+x^2+2y^2}\begin{pmatrix}2x\\\\4y\end{pmatrix}.
$$

At $p=(1,1)$: $1+q(p)=4$, so
$$
\nabla f(p)=\frac14\begin{pmatrix}2\\\\4\end{pmatrix}=\begin{pmatrix}0.5\\\\1\end{pmatrix},
\qquad
\boxed{d_0=\begin{pmatrix}-0.5\\\\-1\end{pmatrix}}.
$$

Also $f(p)=\ln 4\approx 1.3863$.

### Question 2: Fixed step $\alpha=1$

The update is $x_{k+1}=x_k-\alpha\,\dfrac{1}{1+x_k^2+2y_k^2}\begin{pmatrix}2x_k\\\\4y_k\end{pmatrix}$ (it is **not** linear in $x_k$).

**First iteration.**
$$
x_1=\begin{pmatrix}1-0.5\\\\1-1\end{pmatrix}=\boxed{\begin{pmatrix}0.5\\\\0\end{pmatrix}},
\qquad f(x_1)=\ln(1.25)\approx0.2231.
$$

**Second iteration.** $\nabla f(x_1)=\dfrac{1}{1.25}\begin{pmatrix}1\\\\0\end{pmatrix}=\begin{pmatrix}0.8\\\\0\end{pmatrix}$, so
$$
x_2=\begin{pmatrix}0.5-0.8\\\\0\end{pmatrix}=\boxed{\begin{pmatrix}-0.3\\\\0\end{pmatrix}},
\qquad f(x_2)=\ln(1.09)\approx0.0862.
$$

The objective decreases: $1.3863\to0.2231\to0.0862$.

### Question 3: Fixed step $\alpha=2$

$$
x_1=\begin{pmatrix}1-1\\\\1-2\end{pmatrix}=\begin{pmatrix}0\\\\-1\end{pmatrix},\qquad f(x_1)=\ln3\approx1.0986.
$$

Then $\nabla f(x_1)=\dfrac13\begin{pmatrix}0\\\\-4\end{pmatrix}$, so
$$
x_2=\begin{pmatrix}0\\\\-1+\tfrac83\end{pmatrix}=\begin{pmatrix}0\\\\5/3\end{pmatrix},
\qquad f(x_2)=\ln\!\left(1+\tfrac{50}{9}\right)=\ln\tfrac{59}{9}\approx1.8803.
$$

**Observation:** $f$ goes $1.3863\to1.0986\to1.8803$. The first step decreases $f$, the second increases it above $f(x_0)$. The direction is a descent direction at each point, but the step is too long: the iterates jump from one side of the valley to the other.

### Question 4: Exact line search from $p$

$$
x(\alpha)=p+\alpha d_0=\begin{pmatrix}1-0.5\alpha\\\\1-\alpha\end{pmatrix}.
$$

$$
q(\alpha)=(1-0.5\alpha)^2+2(1-\alpha)^2=3-5\alpha+2.25\alpha^2,
$$
so
$$
\boxed{\phi(\alpha)=\ln\!\left(4-5\alpha+2.25\alpha^2\right)},\qquad
\phi'(\alpha)=\frac{-5+4.5\alpha}{4-5\alpha+2.25\alpha^2}.
$$

**Check.** $\phi'(0)=-5/4=-1.25$, and $\|\nabla f(p)\|^2=0.25+1=1.25$, so $\phi'(0)=-\|\nabla f(p)\|^2$ ✓.

The denominator is always positive (discriminant $25-4\cdot4\cdot2.25=-11<0$), so $\phi'(\alpha)=0$ only when the numerator vanishes:
$$
-5+4.5\alpha^*=0\ \Longrightarrow\ \boxed{\alpha_0^*=\tfrac{10}{9}\approx1.111}.
$$

The numerator changes sign from $-$ to $+$, so this is the minimum of $\phi$ along the line.

$$
x_1=\begin{pmatrix}1-\tfrac59\\\\[2pt]1-\tfrac{10}{9}\end{pmatrix}=\boxed{\begin{pmatrix}4/9\\\\-1/9\end{pmatrix}},\qquad
1+q(x_1)=1+\tfrac{16+2}{81}=\tfrac{11}{9},\qquad
f(x_1)=\ln\tfrac{11}{9}\approx0.2007.
$$

(Cross-check: $4-5\cdot\tfrac{10}{9}+2.25\cdot\tfrac{100}{81}=\tfrac{11}{9}$ ✓.)

### Question 5: Orthogonality

$$
\nabla f(x_1)=\frac{9}{11}\begin{pmatrix}8/9\\\\-4/9\end{pmatrix}=\begin{pmatrix}8/11\\\\-4/11\end{pmatrix},
\qquad
\nabla f(x_1)^Td_0=\frac{8}{11}(-0.5)+\left(-\frac{4}{11}\right)(-1)=-\frac4{11}+\frac4{11}=\boxed{0}.
$$

Since $\phi'(\alpha)=\nabla f(x(\alpha))^Td_0$, the condition $\phi'(\alpha^*)=0$ says exactly that **the new gradient is orthogonal to the direction just travelled**: the line is tangent to the level curve through $x_1$. Hence $d_1=-\nabla f(x_1)=\left(-\tfrac8{11},\tfrac4{11}\right)$ is orthogonal to $d_0$, and consecutive steepest-descent directions zig-zag.

### Question 6: Second exact line search

Direction $d_1=\left(-\tfrac8{11},\tfrac4{11}\right)$. Along the line,
$$
x(\alpha)=\begin{pmatrix}\tfrac49-\tfrac8{11}\alpha\\\\[2pt]-\tfrac19+\tfrac4{11}\alpha\end{pmatrix},\qquad
q(\alpha)=\frac29-\frac{80}{99}\alpha+\frac{96}{121}\alpha^2 .
$$

Since $\phi_1=\ln(1+q(\alpha))$ and $1+q>0$, the exact step solves $q'(\alpha)=0$:
$$
-\frac{80}{99}+\frac{192}{121}\alpha=0\ \Longrightarrow\ \boxed{\alpha_1^*=\frac{55}{108}\approx0.5093}.
$$

$$
x_2=\begin{pmatrix}\tfrac49-\tfrac8{11}\cdot\tfrac{55}{108}\\\\[2pt]-\tfrac19+\tfrac4{11}\cdot\tfrac{55}{108}\end{pmatrix}
=\boxed{\begin{pmatrix}2/27\\\\2/27\end{pmatrix}},
\qquad
f(x_2)=\ln\tfrac{247}{243}\approx0.0163.
$$

(Check: $1+q(x_2)=1+3\cdot\tfrac{4}{729}=\tfrac{741}{729}=\tfrac{247}{243}$.)

### Question 7: Comparison

| Method | $f(x_0)$ | $f(x_1)$ | $f(x_2)$ |
|---|---|---|---|
| Fixed step $\alpha=1$ | 1.3863 | 0.2231 | 0.0862 |
| Fixed step $\alpha=2$ | 1.3863 | 1.0986 | 1.8803 |
| Exact line search | 1.3863 | 0.2007 | 0.0163 |

Exact line search gives the lowest value at every iteration, and its step changes from $10/9$ to $55/108$. A fixed step cannot adapt: $\alpha=1$ is acceptable at the start and $\alpha=2$ is too long.

### Question 8: Remark

For a pure quadratic $x^2+2y^2$, the update is $(x,y)\mapsto((1-2\alpha)x,(1-4\alpha)y)$, which is linear, and convergence requires $|1-2\alpha|<1$ and $|1-4\alpha|<1$, i.e. $\alpha<1/2$.

Here the gradient is divided by $1+q$, so the update is nonlinear and no single global condition exists: the effective step is $\alpha/(1+q)$, which changes along the path. Near the minimizer $q\approx0$, the factor tends to $1$, $f\approx x^2+2y^2$, and the local condition is again $\alpha<1/2$.

**Warning on $\alpha=1$:** it worked in Question 2 only because $y_1=0$ exactly. Any small perturbation of $y$ near the origin is multiplied by $1-4\alpha=-3$ and grows. So a fixed step that "works" on one trajectory can still be unstable.

---

# Exercise 2: Armijo backtracking when exact line search is not available

Let
$$
f(x,y)=x^4+y^2,\qquad p=(1,1),\qquad d=-\nabla f(p).
$$

**Aim.** Understand Armijo as a rule that accepts or rejects a step by comparing the actual decrease of $f$ with a fraction $c_1$ of the decrease predicted by the tangent line. See the roles of $c_1$ and $\beta$, and see that no equation has to be solved.

## Statement

1. Compute $\nabla f(p)$, $d$, then $\phi(\alpha)=f(p+\alpha d)$, $\phi(0)$ and $\phi'(0)$. Check that $\phi'(0)=-\|\nabla f(p)\|^2$.
2. Write the equation $\phi'(\alpha)=0$. Why can it not be solved by hand? Locate its root $\alpha^*$ to three decimals (sign table), and give $\phi(\alpha^*)$.
3. Write the Armijo condition for $c_1=0.2$. Run backtracking from $\alpha_0=1$ with $\beta=0.5$. For each trial, give $\phi(\alpha)$, the threshold, and the decision.
4. Give the accepted point $x_1$ and $f(x_1)$. Compare with exact line search: what fraction of the best possible decrease along the line is obtained?
5. Repeat with $c_1=0.01$ and $c_1=0.9$ (same $\alpha_0$ and $\beta$). Compare the number of trials, the accepted step and the decrease obtained. Interpret.
6. *(Bonus)* Find the largest $c_1$ for which the exact step $\alpha^*$ is accepted.

## Solution

### Question 1: Direction and $\phi$

$$
\nabla f(x,y)=\begin{pmatrix}4x^3\\\\2y\end{pmatrix},\qquad
\nabla f(p)=\begin{pmatrix}4\\\\2\end{pmatrix},\qquad
\boxed{d=\begin{pmatrix}-4\\\\-2\end{pmatrix}}.
$$

$$
x(\alpha)=\begin{pmatrix}1-4\alpha\\\\1-2\alpha\end{pmatrix},\qquad
\boxed{\phi(\alpha)=(1-4\alpha)^4+(1-2\alpha)^2}.
$$

$$
\phi'(\alpha)=-16(1-4\alpha)^3-4(1-2\alpha),\qquad
\phi(0)=2,\qquad \phi'(0)=-16-4=-20 .
$$

Check: $\|\nabla f(p)\|^2=16+4=20$, so $\phi'(0)=-\|\nabla f(p)\|^2$ ✓.

### Question 2: Exact step is not available by hand

$\phi'(\alpha)=0$ is equivalent to $4(1-4\alpha)^3+(1-2\alpha)=0$, i.e.
$$
256\alpha^3-192\alpha^2+50\alpha-5=0,
$$
a cubic with no simple rational root. (Here $\phi''(\alpha)=192(1-4\alpha)^2+8>0$, so $\phi$ is convex and the root is unique, but finding it needs a numerical method.)

Sign table:

| $\alpha$ | 0.35 | 0.354 | 0.355 | 0.36 |
|---|---|---|---|---|
| $\phi'(\alpha)$ | $-0.176$ | $-0.016$ | $+0.025$ | $+0.243$ |

So $\boxed{\alpha^*\approx0.354}$ and $\phi(\alpha^*)\approx0.115$. The best possible decrease along the line is $2-0.115=1.885$.

### Question 3: Armijo with $c_1=0.2$, $\beta=0.5$

The condition is $\phi(\alpha)\le\phi(0)+c_1\alpha\phi'(0)$, i.e.
$$
\phi(\alpha)\le 2-0.2\cdot20\,\alpha=2-4\alpha .
$$

| Trial $\alpha$ | $\phi(\alpha)$ | Threshold $2-4\alpha$ | Decision |
|---|---|---|---|
| 1 | $81+1=82$ | $-2$ | Rejected |
| 0.5 | $1+0=1$ | $0$ | Rejected |
| 0.25 | $0+0.25=0.25$ | $1$ | **Accepted** |

Note that at $\alpha=0.5$ the point is $(-1,0)$: it has passed over the valley in the $x$ direction and $f=1$ is not below the threshold $0$.

### Question 4: Accepted point and comparison

$$
x_1=p+0.25\,d=\boxed{\begin{pmatrix}0\\\\0.5\end{pmatrix}},\qquad f(x_1)=0.25 .
$$

Exact search gives $\phi(\alpha^*)\approx0.115$. The Armijo step decreases $f$ by $2-0.25=1.75$, i.e. about $\dfrac{1.75}{1.885}\approx93\%$ of the best possible decrease along the line, using three evaluations of $\phi$ and no root-finding.

### Question 5: Effect of $c_1$

| $c_1$ | Threshold | Trials ($\alpha$) | Accepted $\alpha$ | $f(x_1)$ | Decrease (share of 1.885) | Accepted set |
|---|---|---|---|---|---|---|
| 0.01 | $2-0.2\alpha$ | 1, 0.5 (2 trials) | 0.5 | 1 | 1.00 (53%) | $[0,\,0.543]$ |
| 0.2 | $2-4\alpha$ | 1, 0.5, 0.25 (3 trials) | 0.25 | 0.25 | 1.75 (93%) | $[0,\,0.429]$ |
| 0.9 | $2-18\alpha$ | 1, 0.5, …, 0.0156 (7 trials) | 0.0156 | 1.711 | 0.289 (15%) | $[0,\,0.021]$ |

**Interpretation.**

- $c_1=0.01$: the condition is very permissive. It accepts $\alpha=0.5$, which has already overshot the minimizer $\alpha^*\approx0.354$ and gives only about half of the possible decrease. Armijo guards against bad steps, not against mediocre ones.
- $c_1=0.9$: the condition demands almost the full decrease predicted by the tangent line. Since $\phi$ is strongly convex, it rises above its tangent quickly, so only tiny steps pass. Backtracking needs 7 trials and ends with a very short step that excludes $\alpha^*$.
- $c_1=0.2$ is a balanced case here. In practice a small value such as $10^{-4}$ is used, and an additional condition (the Wolfe curvature condition) is added to forbid steps that are too short.
- $\beta$ controls how fast the step shrinks: a smaller $\beta$ needs fewer trials but may jump over good steps.

### Question 6 (bonus): Largest $c_1$ accepting $\alpha^*$

$\alpha^*$ is accepted if $\phi(\alpha^*)\le 2-20c_1\alpha^*$, i.e.
$$
c_1\le\frac{2-\phi(\alpha^*)}{20\,\alpha^*}\approx\frac{1.885}{7.08}\approx0.266 .
$$

For any $c_1>0.266$ the exact minimizer along the line is rejected, which shows that Armijo does not look for the minimum: it only defines an acceptable range of steps.

---

## Summary of key numerical results

| Quantity | Exercise 1 | Exercise 2 |
|---|---|---|
| $f(p)$ | $\ln4\approx1.386$ | $2$ |
| $\nabla f(p)$ | $(0.5,\ 1)$ | $(4,\ 2)$ |
| $\phi'(0)$ | $-1.25$ | $-20$ |
| Exact step | $10/9$ then $55/108$ | $\approx0.354$ (cubic) |
| Armijo step ($c_1=0.2$, $\beta=0.5$) | not used | $0.25$ |
