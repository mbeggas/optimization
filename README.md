# A/ Convexity, Taylor Approximation and Line Search
## 1. Convexity

Convexity determines whether an optimization problem is computationally tractable. In non-convex optimization, algorithms can get trapped in local minima, saddle points, or flat regions. In convex optimization, **any local minimum is guaranteed to be a global minimum**.

### Mathematical Definition

A set $C \subseteq \mathbb{R}^d$ is convex if the line segment connecting any two points in $C$ lies entirely inside $C$.

A function $f: C \to \mathbb{R}$ defined on a convex set $C$ is convex if, for all $x, y \in C$ and $\lambda \in [0, 1]$:

$$f(\lambda x + (1 - \lambda)y) \le \lambda f(x) + (1 - \lambda) f(y)$$

Geometrically, this means the line segment (secant line) connecting $(x, f(x))$ and $(y, f(y))$ lies on or above the graph of $f$.

### Characterizations for Differentiable Functions

* **First-Order Condition:** A continuously differentiable function $f$ is convex if and only if for all $x, y$:

$$f(y) \ge f(x) + \nabla f(x)^T (y - x)$$



This means the first-order Taylor approximation at any point $x$ serves as a **global lower bound** for the function.
* **Second-Order Condition:** A twice continuously differentiable function $f$ is convex if and only if its Hessian matrix $\nabla^2 f(x)$ is positive semi-definite ($\nabla^2 f(x) \succeq 0$) for all $x$:

$$v^T \nabla^2 f(x) v \ge 0 \quad \forall v \in \mathbb{R}^d$$



This indicates that the function has non-negative curvature in every direction.

---

## 2. Taylor Approximation

Optimization algorithms locally approximate complex objective functions using polynomials via **Taylor expansions**. These local approximations govern how updates and search directions are chosen.

Given a twice-differentiable function $f: \mathbb{R}^d \to \mathbb{R}$, its multivariate Taylor series expansion around point $x_k$ for a perturbation vector $\Delta x$ is:

$$f(x_k + \Delta x) = f(x_k) + \nabla f(x_k)^T \Delta x + \frac{1}{2} \Delta x^T \nabla^2 f(x_k) \Delta x + \mathcal{O}(\Vert{}\Delta x\Vert{}^3)$$

### First-Order Taylor Approximation (Linear Model)

Truncating the series after the linear term gives:

$$f(x_k + \Delta x) \approx f(x_k) + \nabla f(x_k)^T \Delta x$$

* **Application in Gradient Descent:** To minimize this local linear model for a step of fixed length $\Vert{}\Delta x\Vert{} = \alpha$, the step direction must minimize the inner product $\nabla f(x_k)^T \Delta x$. By the Cauchy-Schwarz inequality, this occurs when $\Delta x = -\alpha \frac{\nabla f(x_k)}{\Vert{}\nabla f(x_k)\Vert{}}$, yielding the standard gradient descent step.

### Second-Order Taylor Approximation (Quadratic Model)

Including the quadratic term yields:

$$f(x_k + \Delta x) \approx f(x_k) + \nabla f(x_k)^T \Delta x + \frac{1}{2} \Delta x^T \nabla^2 f(x_k) \Delta x$$

* **Application in Newton's Method:** Newton's method finds the optimal step $\Delta x$ by directly setting the gradient of this quadratic model with respect to $\Delta x$ to zero:

$$\nabla f(x_k) + \nabla^2 f(x_k) \Delta x = 0 \implies \Delta x = -[\nabla^2 f(x_k)]^{-1} \nabla f(x_k)$$



Because it accounts for local curvature via the Hessian $\nabla^2 f(x_k)$, Newton's method achieves quadratic convergence near a local minimum.

---

## 3. Line Search

Given a descent direction $p_k$ (such as $p_k = -\nabla f(x_k)$ in gradient descent or $p_k = -[\nabla^2 f(x_k)]^{-1}\nabla f(x_k)$ in Newton's method), **line search** is the procedure used to select an optimal step size $\alpha_k > 0$ for updating $x_{k+1} = x_k + \alpha_k p_k$.

Define the 1D objective function along the search ray as:

$$\phi(\alpha) = f(x_k + \alpha p_k)$$

### Exact Line Search

Exact line search computes the step size that globally minimizes $f$ along the ray:

$$\alpha_k = \arg\min_{\alpha > 0} f(x_k + \alpha p_k)$$

* **Limitation:** Solving a 1D optimization problem exactly at every iteration is computationally prohibitive for expensive high-dimensional objective functions.

### Inexact Line Search & The Armijo Condition

In practice, algorithms search for a step size $\alpha_k$ that guarantees a **sufficient decrease** in $f$ without taking steps that are too small.

The **Armijo Condition** enforces that the reduction in $f$ is proportional to both step size $\alpha$ and the local directional derivative $\nabla f(x_k)^T p_k$:

$$f(x_k + \alpha p_k) \le f(x_k) + c_1 \alpha \nabla f(x_k)^T p_k \quad \text{where } c_1 \in (0, 1)$$

Because $p_k$ is a descent direction, $\nabla f(x_k)^T p_k < 0$, making $c_1 \alpha \nabla f(x_k)^T p_k$ a negative offset from $f(x_k)$.

### Backtracking Line Search Algorithm

Backtracking provides an efficient iterative approach to find a valid step size $\alpha$:

1. **Parameters:** Choose initial step $\alpha = 1$, contraction factor $\rho \in (0, 1)$ (e.g., $\rho = 0.5$), and constant $c_1 \in (0, 1)$ (e.g., $c_1 = 10^{-4}$).
2. **Evaluation:** Check if $f(x_k + \alpha p_k) \le f(x_k) + c_1 \alpha \nabla f(x_k)^T p_k$.
3. **Update Step:** If the condition fails, reduce step size:

$$\alpha \leftarrow \rho \cdot \alpha$$



Repeat step 2 until the condition holds.
4. **Set Step:** Assign $\alpha_k = \alpha$ and update $x_{k+1} = x_k + \alpha_k p_k$.

-----

# B/ Gradient Descent

## 1. Problem Formulation & Core Intuition

Gradient Descent is a first-order iterative optimization algorithm designed to solve unconstrained optimization problems of the form:

$$\min_{x \in \mathbb{R}^d} f(x)$$

where $f: \mathbb{R}^d \to \mathbb{R}$ is continuously differentiable.

The core objective is to generate a sequence of iterates $x_0, x_1, x_2, \dots \in \mathbb{R}^d$ such that $f(x_{k+1}) < f(x_k)$, eventually converging to a point $x^*$ where the gradient vanishes: $\nabla f(x^*) = 0$.

---

## 2. Derivation via First-Order Taylor Expansion

The algorithm's update rule is derived directly from the **first-order Taylor approximation** of $f(x)$.

Suppose we are at iterate $x_k$ and want to take a small step $\Delta x$ of length $\Vert{}\Delta x\Vert{} = \alpha > 0$. Expanding $f(x_k + \Delta x)$ locally:

$$f(x_k + \Delta x) \approx f(x_k) + \nabla f(x_k)^T \Delta x$$

To achieve the maximum decrease in $f(x)$, we want to choose $\Delta x$ to make the inner product term $\nabla f(x_k)^T \Delta x$ as negative as possible.

Let $\Delta x = \alpha d$, where $d \in \mathbb{R}^d$ is a unit search direction ($\Vert{}d\Vert{} = 1$). Expanding the inner product using the definition of the dot product yields:

$$\nabla f(x_k)^T \Delta x = \alpha \Vert{}\nabla f(x_k)\Vert{} \Vert{}d\Vert{} \cos(\theta) = \alpha \Vert{}\nabla f(x_k)\Vert{} \cos(\theta)$$

where $\theta$ is the angle between the gradient vector $\nabla f(x_k)$ and the step direction $d$.

* The value of $\cos(\theta)$ is minimized when $\theta = 180^\circ$ ($\cos(\theta) = -1$).
* Consequently, the direction $d$ that minimizes the local approximation is directly opposite to the gradient:

$$d = -\frac{\nabla f(x_k)}{\Vert{}\nabla f(x_k)\Vert{}}$$



Multiplying direction $d$ by step size $\alpha$ gives the **Gradient Descent Update Rule**:

$$x_{k+1} = x_k - \alpha_k \nabla f(x_k)$$

---

## 3. Convergence Analysis & The Role of Convexity

The behavior and speed of gradient descent depend heavily on whether $f(x)$ is **convex** and possesses a **Lipschitz continuous gradient**.

### $L$-Smoothness & The Descent Lemma

A function $f$ has an $L$-Lipschitz continuous gradient ($L$-smooth) if for all $x, y$:

$$\Vert{}\nabla f(x) - \nabla f(y)\Vert{} \le L \Vert{}x - y\Vert{}$$

This places an upper bound $L$ on the function's curvature. From this, we obtain the quadratic upper bound known as the **Descent Lemma**:

$$f(y) \le f(x) + \nabla f(x)^T(y - x) + \frac{L}{2} \Vert{}y - x\Vert{}^2$$

Substituting the gradient descent step $y = x_k - \alpha \nabla f(x_k)$:

$$f(x_{k+1}) \le f(x_k) - \alpha \Vert{}\nabla f(x_k)\Vert{}^2 + \frac{L \alpha^2}{2} \Vert{}\nabla f(x_k)\Vert{}^2 = f(x_k) - \alpha \left(1 - \frac{\alpha L}{2}\right) \Vert{}\nabla f(x_k)\Vert{}^2$$

To guarantee that $f(x_{k+1}) < f(x_k)$, we require $\alpha \left(1 - \frac{\alpha L}{2}\right) > 0$, which holds when $0 < \alpha < \frac{2}{L}$.

Maximizing this decrease gives the optimal fixed step size $\alpha = \frac{1}{L}$, yielding:

$$f(x_{k+1}) \le f(x_k) - \frac{1}{2L} \Vert{}\nabla f(x_k)\Vert{}^2$$

### Convergence Rates by Function Class

| Function Class | Curvature Condition | Convergence Rate | Iterations to reach $\epsilon$-accuracy |
| --- | --- | --- | --- |
| **Convex ($L$-Smooth)** | $f(y) \ge f(x) + \nabla f(x)^T(y-x)$ | $\mathcal{O}(1/k)$ (Sublinear) | $\mathcal{O}(1/\epsilon)$ |
| **Strongly Convex ($\mu, L$)** | $\mu I \preceq \nabla^2 f(x) \preceq L I$ | $\mathcal{O}\left(\left(1 - \frac{\mu}{L}\right)^k\right)$ (Linear) | $\mathcal{O}\left(\kappa \log(1/\epsilon)\right)$ |
| **Non-Convex ($L$-Smooth)** | Arbitrary curvature | Stationary Point: $\min_{i \le k} \Vert{}\nabla f(x_i)\Vert{} \le \mathcal{O}(1/\sqrt{k})$ | $\mathcal{O}(1/\epsilon^2)$ |

*Note on Condition Number:* For strongly convex functions, $\kappa = \frac{L}{\mu} \ge 1$ is the condition number. High condition numbers ($\kappa \gg 1$) cause gradient descent to oscillate heavily across steep valleys rather than moving down the center.

---

## 4. Step Size Selection via Line Search

Because the Lipschitz constant $L$ is rarely known in practice, static step sizes can cause divergence or slow progress. **Inexact Line Search** solves this by dynamically adapting $\alpha_k$ at every step.

Substituting the gradient descent search direction $p_k = -\nabla f(x_k)$ into the **Armijo Condition**:

$$f(x_k - \alpha \nabla f(x_k)) \le f(x_k) - c_1 \alpha \Vert{}\nabla f(x_k)\Vert{}^2 \quad \text{where } c_1 \in (0, 1)$$

### Integrated Backtracking Algorithm

1. Start with initial step $\alpha = 1$, parameters $\rho = 0.5$, $c_1 = 10^{-4}$.
2. Evaluate:

$$\text{Check if } f(x_k - \alpha \nabla f(x_k)) \le f(x_k) - c_1 \alpha \Vert{}\nabla f(x_k)\Vert{}^2$$


3. If condition fails, shrink step: $\alpha \leftarrow \rho \cdot \alpha$.
4. Update iterate: $x_{k+1} = x_k - \alpha \nabla f(x_k)$.

This automatically finds an adaptive step size that respects the local upper curvature bound without needing global information about $L$.

---

## 5. Main Variants of Gradient Descent

In machine learning and empirical risk minimization, the objective function is typically a sum over dataset samples $N$:

$$f(x) = \frac{1}{N} \sum_{i=1}^{N} f_i(x)$$

### Batch Gradient Descent (BGD)

Computes the exact gradient over the entire dataset:


$$x_{k+1} = x_k - \alpha \left(\frac{1}{N} \sum_{i=1}^{N} \nabla f_i(x_k)\right)$$

* **Pros:** Deterministic convergence guarantees.
* **Cons:** Extremely expensive when $N$ is large ($O(N)$ per iteration).

### Stochastic Gradient Descent (SGD)

Approximates the gradient using a single randomly chosen sample index $i_k \in \{1, \dots, N\}$:


$$x_{k+1} = x_k - \alpha \nabla f_{i_k}(x_k)$$

* **Pros:** $O(1)$ computation per iteration; noise helps escape saddle points in non-convex landscapes.
* **Cons:** High variance requires decaying step sizes $\alpha_k \to 0$ to ensure convergence.

### Mini-Batch Gradient Descent

Computes the gradient over a random subset $B \subset \{1, \dots, N\}$ of size $\vert{}B\vert{} = b$:


$$x_{k+1} = x_k - \alpha \left(\frac{1}{b} \sum_{i \in B} \nabla f_i(x_k)\right)$$

* Balances the computational speed of SGD with the vectorization efficiency and reduced variance of Batch GD.
