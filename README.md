# uncoupled
R package implementing minimum-Wasserstein isotonic regression estimator from Rigollet and Weed (2019)

# Implementation details

To implement [Rigollet and Weed](https://arxiv.org/abs/1806.10648)'s relaxed minimum-$W_2$ estimator, we find

$$
\begin{aligned}
\hat{\mu} \in \text{argmin}_{\mu \in \mathcal{M}_{\mathcal{A}, V}} W_2^2 \left(\Pi_{\mathcal{A}\sharp}(\mu \star \mathcal{D}), \hat{\pi} \right)
\end{aligned}
$$

where $\Pi_{\mathcal{A}\sharp}(\mu \star \mathcal{D})$ is the pushforward of the measure $\mu \star \mathcal{D}$ by the projection operator defined in [Rigollet and Weed](https://arxiv.org/abs/1806.10648). Because both pushforward measures considered here are discrete, finding the Wasserstein distance between a given $\Pi_{\mathcal{A}\sharp}(\mu \star \mathcal{D})$ and $\hat{\pi}$ is implemented as a linear program

$$
\begin{aligned}
\text{minimize}~~\sum_{ij} \gamma_{ij} &(y_i - \alpha_j)^2 \\
\text{subject to}~~\sum_{i} \gamma_{ij} &= a_j &j = 0, \dots, N \\
\sum_{j} \gamma_{ij} &= \hat{\pi}_i &i = 1, \dots, n \\
\gamma_{ij} &\geq 0
\end{aligned}
$$

where $\alpha_j$ is the $j^{\text{th}}$ partition of the real line defined in [Rigollet and Weed](https://arxiv.org/abs/1806.10648), $y_i$ is the $i^{\text{th}}$ observed data point, and $a_j = \Pi_{\mathcal{A}\sharp}(\mu \star \mathcal{D})_j$.

From the dual problem to this linear program, we obtain the subgradient of the $W_2^2$ distance with respect to the measure $a$. This subgradient is equal to the vector of optimal dual solutions corresponding to the first equality constraint in the optimization problem above. So we have

$$
\begin{aligned}
\frac{\partial W_2^2}{\partial a_j} = g_j^*~~~~~~~~j = 0, \dots, N \\
\end{aligned}
$$

where $g_j^*$ is the $j^{\text{th}}$ component of the optimal dual vector. Next, define $\mu_i := \mu (\alpha_i)$; in other words the measure weight given to $\alpha_i$. Then based on the definition of the projection operator $\Pi_{\mathcal{A}\sharp}$, we have

$$
\begin{aligned}
a_j (\mu_1, \dots, \mu_N) = \sum_{i=1}^N \mu_i \int_{\alpha_j}^{\alpha_{j+1}} p(x ; \alpha_j, \sigma^2) dx
\end{aligned}
$$

where the noise distribution $\mathcal{D}$ has density function $p$. We compute partial derivatives of the $a$'s with respect to each $\mu_i$ to obtain

$$
\begin{aligned}
\frac{\partial a_j}{\partial \mu_i} = \int_{\alpha_j}^{\alpha_{j+1}} p(x ; \alpha_j, \sigma^2) dx
\end{aligned}
$$

Next we obtain subgradients of the $W_2^2$ objective with respect to the measure weights $\mu_1, \dots, \mu_N$ with

$$
\begin{aligned}
\frac{\partial W_2^2}{\partial \mu_i} = \sum_{j} \frac{\partial W_2^2}{\partial a_j} \frac{\partial a_j}{\partial \mu_i}
\end{aligned}
$$

and use projected subgradient updates of the form

$$
\begin{aligned}
\mu_{t+1} = \text{Proj}_{\mu \in \mathcal{M}_{\mathcal{A}, V}} \left( \mu_t - \lambda_t \frac{d W_2^2}{d \mu} \right)
\end{aligned}
$$

for step size $\lambda_t$ to arrive at the optimal value of $\mu$. To speed convergence to the optimal $\mu^*$, we use the CFM step-size rule proposed by Camerini P.M., Fratta L., Maffioli F. (1975). The operator $\text{Proj}_{\mu \in \mathcal{M}_{\mathcal{A}, V}}$ is the projection operator onto the probability simplex.
