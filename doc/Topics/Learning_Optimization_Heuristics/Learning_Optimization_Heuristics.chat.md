https://chatgpt.com/share/6af042f7-1e18-83eb-9768-61a3854ce1f1

---

# USER

For long time I was interested in making effiicent optimization algorithms resp heuristics but simple, not like BFGS anr CG which are mathematicalyl sound and more or less optimal close to parabolic minimum, but something which can effciently walk complex non-linear no convex surface (long curved valleays and canoyons) using simple herustic using either only energy or energy ang gradient (force). Lets discus them separately but the ideas are simple.

Typically these algorithms use some kide of intertia or memoy mixing old gradient. typical example is FIRE algorithm (fast interial relaxation engine) [https://arxiv.org/pdf/1908.02038](https://arxiv.org/pdf/1908.02038)
or ADAM used for machine lerning ([https://arxiv.org/abs/1412.6980](https://arxiv.org/abs/1412.6980)), then we have gradient less optimizations like POWELL method&#x20;
[https://en.wikipedia.org/wiki/Powell%27s_method](https://en.wikipedia.org/wiki/Powell%27s_method)
[https://en.wikipedia.org/wiki/Successive_parabolic_interpolation](https://en.wikipedia.org/wiki/Successive_parabolic_interpolation)

What I was using often is something which is like gradient-gree ADAM where step is linear combination of random exploration previous successful step (which decay at every unsuccessful atempt), or that we iteratively build some low-dimensional aporximation of gradient or hessian by every unsucsefull move...&#x20;

Now lets approch it as machine learning heuristics rather than rigorous mathematical problem with rigorous asumption. We have a generic potential energy surface no which we want to walk toward minimum, but we do not know what it is. We want to have some simple heuristic/machine-lerning algorithm which learns small number of parameters how to walk on it. The algorithm should be fast and simple (if we use it e.g. for classical forcefield, the machinery of the optimization algorithm shold be computationally much cheapr that the evaluation of the energy and force using classical forcediel). Not sure we should use neural-network/multi-layer-preceptron with just few neurons/layers (e.g. <10 layes each with <10 neurons ?), or some symbolic regression ...&#x20;

Also question is if we should do learning on the fly, or rather pre train optimal parameters on large number of problems  from similar class.&#x20;

Ideal would be if the network is pre-tranined on many similar problems but then it can quickly fine tune itself for given potential using this pretrained prior...

Are there any works like this? I guess there should be. But really preferentially somethin simple which does not use much computing resources during uses (it can use large resources to pretran, e.g. terabytes of molecular dynamics trajectroies with forces in each step, but not much neurons during inference when doing the actual minimization steps)

---
---

# Chat GPT 6


Yes. There is a substantial research field called Learning to Optimize (L2O) or meta-optimization, and several papers are remarkably close to what you propose. However, much of this research uses unnecessarily complicated neural networks for your particular goal.

I think the most interesting direction is not to train a neural network to predict the next atomic configuration, but to learn a very small optimization algorithm that uses the history of forces, displacements, and successful or unsuccessful steps to navigate an unknown potential-energy surface.

The distinction matters: the network should learn how to optimize, not memorize the potential itself.

## 1. Existing research closest to your idea

These are the papers I would start with, ranked by conceptual relevance rather than citation count.

[Learning to Learn without Gradient Descent by Gradient Descent](https://proceedings.mlr.press/v70/chen17e.html?utm_source=chatgpt.com)

Chen et al., ICML 2017 — gradient-free optimization

A recurrent neural network learns how to propose new evaluation points using the history of previous queries and objective values. The authors train on simple synthetic functions and demonstrate transfer to more complicated black-box optimization problems.

This is particularly close to your gradient-free ADAM with memory of successful moves. The difference is that you want a much simpler update rule, preferably with computational cost \\(O(N)\\).

[Learning to Optimize Molecular Geometries Using Reinforcement Learning](https://pubmed.ncbi.nlm.nih.gov/33470813/?utm_source=chatgpt.com)

Journal of Chemical Theory and Computation, 2021 — molecular geometry optimization

Directly addresses molecular relaxation. A reinforcement-learning model learns corrections to BFGS steps, targeting situations where the quadratic approximation is poor. The authors report approximately 30% fewer optimization steps than conventional BFGS on perturbed molecular geometries.

It is one of the clearest demonstrations that learned update rules can improve genuine molecular optimization.

[Learning to Learn by Gradient Descent by Gradient Descent](https://arxiv.org/abs/1606.04474?utm_source=chatgpt.com)

Andrychowicz et al., NeurIPS 2016 — gradient-based optimization

The classic learned optimizer: replace a handcrafted gradient update with an LSTM that remembers previous gradients. Importantly, the same small network operates independently on every coordinate, making its parameters independent of the problem dimension.

It establishes the basic architecture, although a coordinate-wise optimizer has limitations for strongly coupled coordinates.

[Symbolic Discovery of Optimization Algorithms](https://arxiv.org/abs/2302.06675?utm_source=chatgpt.com)

Chen et al., NeurIPS 2023 — discovering simple rules

Instead of deploying a neural network, the researchers searched over mathematical programs and discovered Lion, an extremely simple momentum-based optimizer.

This is an excellent example of your desired philosophy: spend substantial resources during algorithm discovery, but produce something consisting of a few arithmetic operations that can be implemented efficiently.

The first and fourth papers are particularly interesting together. One learns a stateful optimization policy, while the other shows that the result of algorithm discovery can be a simple mathematical expression rather than a neural network.

## 2. What exactly should the learned algorithm learn?

I would distinguish three different kinds of memory and adaptation.

| Level                    | What is learned                                                        | When                                 | Computational cost during optimization              |
| ------------------------ | ---------------------------------------------------------------------- | ------------------------------------ | --------------------------------------------------- |
| Offline meta-learning    | General optimization rules and their parameters                        | Once, across many potential surfaces | Essentially zero beyond evaluating the learned rule |
| Online state adaptation  | Momentum, characteristic step size, local curvature, useful directions | Every optimization step              | \\(O(N)\\) or less                                  |
| Online weight adaptation | Changing the neural network's weights for the particular potential     | During each optimization             | Potentially significant                             |

Here \\(N\\) is the number of optimization coordinates.

My preference would be to combine the first two, without initially implementing the third.

This actually gives you much of what you want from pretrained priors followed by rapid fine-tuning. A recurrent optimizer can have fixed weights but adaptive internal state. As it explores a particular potential, that state gradually becomes specialized to its geometry.

One can think of it as:

\\[ \boxed{ \begin{aligned} \mathbf x\_{t+1}&=\mathbf x_t+\Delta\mathbf x_t\\\ (\Delta\mathbf x_t,\mathbf h\_{t+1})  &=\mathcal A\_\theta(\mathbf F_t,E_t,\mathbf h_t) \end{aligned} } \\]

where \\(\theta\\) is learned offline, while \\(\mathbf h_t\\) is updated online.

The practical question is how to choose a minimal but informative \\(\mathbf h_t\\).

## 3. Force-based optimization: following curved valleys

Consider a narrow, curved valley.

Conceptual illustration: an optimizer should preserve motion along the valley while correcting its direction as the valley bends. The drawn lines are illustrative, not simulated trajectories.

The conventional gradient does not necessarily point in the direction in which we want to travel. In a stiff valley, a small transverse displacement creates a strong restoring force. Meanwhile, the longitudinal force can be relatively weak.

This creates two competing requirements:

- Retain the useful longitudinal motion rather than following every instantaneous change in force.
- Correct the trajectory sufficiently quickly when the valley curves, without oscillating between its walls.

FIRE solves part of this by mixing the velocity with the force and adapting the timestep using the sign of \\(\mathbf F\cdot\mathbf v\\). FIRE 2.0 also refines the integration and overshoot handling.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://doi.org\&sz=32)

DOI+1



However, the optimizer could potentially make better decisions using a little more information than FIRE retains.

### A minimal learned generalization of FIRE

Suppose we keep the current force \\(\mathbf F_t\\), the previous step \\(\mathbf s_t=\mathbf x_t-\mathbf x\_{t-1}\\), and one momentum vector \\(\mathbf m_t\\).

We can construct a small set of dimensionless features:

\\[ c_t=\frac{\mathbf F_t\cdot\mathbf m_t} {\\|\mathbf F_t\\|\\|\mathbf m_t\\|+\epsilon} \\]

This measures the alignment between the current force and the remembered direction.

Another useful feature is force-direction persistence:

\\[ q_t=\frac{\mathbf F_t\cdot\mathbf F\_{t-1}} {\\|\mathbf F_t\\|\\|\mathbf F\_{t-1}\\|+\epsilon}. \\]

And perhaps most interestingly, estimate the curvature along the previous step:

\\[ k_t \approx -\frac{\mathbf s_t\cdot(\mathbf F_t-\mathbf F\_{t-1})} {\\|\mathbf s_t\\|^2+\epsilon}. \\]

The last expression is effectively a directional secant approximation to the Hessian, obtained without additional force evaluations.

Rather than feed all atomic coordinates to a network, use these few scalars to control an optimizer such as

\\[ \begin{aligned} \mathbf m\_{t+1}&=\beta_t\mathbf m_t+(1-\beta_t)\mathbf F_t,\\\ \Delta\mathbf x_t&=\eta_t\mathbf m\_{t+1}. \end{aligned} \\]

The learned controller predicts

\\[ (\beta_t,\eta_t)= \operatorname{MLP}\_\theta(c_t,q_t,k_t,\ldots). \\]

For example, a network with six inputs, eight hidden neurons and two outputs has only 74 trainable weights and biases.

More importantly, the neural network is evaluated only once per optimization step, not once per atom. It chooses a few scalar parameters controlling vector operations.

Such a controller could learn to increase momentum in persistent, gently curving valleys and reduce it when the force reverses or the local curvature changes rapidly.

This is a proposed architecture, not a published algorithm whose superior performance has been demonstrated.

### Why the scalar-controller architecture matters

A generic coordinate-wise MLP, like a small ADAM replacement, cannot directly understand that two coordinate components belong to the same stiff collective mode. It processes them independently.

A controller operating on dot products, norms and history vectors instead respects rotations of Cartesian coordinates. It can also be applied to systems with different numbers of atoms.

There is a trade-off: a few global scalars cannot adapt separate regions of a molecule with very different stiffness. One possible extension is to run the same controller independently for different blocks: rigid-body translation, rotation, internal molecular coordinates, substrate degrees of freedom, etc.

For your intended small classical-force-field systems, I would test the global version first.

## 4. Energy-only optimization: learning from successful and unsuccessful moves

I actually find your gradient-free idea especially interesting, because it allows an extremely lightweight implementation.

There is also a very relevant precedent that is easy to overlook: [Covariance Matrix Adaptation Evolution Strategies (CMA-ES)](https://cma-es.github.io/?utm_source=chatgpt.com).

CMA-ES remembers successful directions, accumulates an evolution path, and gradually changes the distribution of proposed displacements. The [(1+1)-CMA-ES](https://www.cmap.polytechnique.fr/~nikolaus.hansen/ACECMUaa1p1CMAfES.pdf?utm_source=chatgpt.com) version is particularly close to your idea of maintaining a single current configuration, generating trial moves and adapting from their success or failure.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://cma-es.github.io\&sz=32)

cma-es.github.io+1



However, full covariance adaptation generally has \\(O(N^2)\\) memory and arithmetic overhead. You probably don't want that.

### A. A two-vector optimizer

Let us try an extremely simple approach.

We maintain just two vectors:

- \\(\mathbf p\\): accumulated successful displacement, analogous to momentum.
- \\(\boldsymbol\xi\\): a freshly sampled random exploration direction.

Propose the trial displacement

\\[ \boxed{\Delta\mathbf x =a\mathbf p+b\boldsymbol\xi} \\]

where \\(a,b\\) are adaptive scalars.

After evaluating the trial energy:

\\[ \mathbf p\_{t+1}= \begin{cases} \beta\mathbf p_t+(1-\beta)\Delta\mathbf x_t, & E\_{\rm trial}\<E_t,\\\ \gamma\mathbf p_t, & E\_{\rm trial}\ge E_t. \end{cases} \\]

Here \\(\gamma\\) is a decay factor on rejected steps.

This is almost exactly the family of heuristics you described.

Even without learning, it is a plausible optimizer. But it has a serious limitation: a rejected move only tells us to reduce the remembered direction. We throw away information about which component of the proposed displacement was wrong.

### B. Treat every rejected move as a measurement of the potential

Suppose we propose a displacement \\(\mathbf s\\) and observe

\\[ \Delta E=E(\mathbf x+\mathbf s)-E(\mathbf x). \\]

For sufficiently small \\(\mathbf s\\),

\\[ \Delta E\approx \mathbf g\cdot\mathbf s, \qquad \mathbf g=\nabla E. \\]

So every evaluation gives us one scalar projection of the unknown gradient.

This works whether the move was accepted or rejected.

We can maintain an estimated gradient and update it using a normalized prediction error:

\\[ \boxed{ \hat{\mathbf g}\_{t+1}= \hat{\mathbf g}\_{t} +\lambda \frac{\Delta E-\hat{\mathbf g}\_{t}\cdot\mathbf s} {\\|\mathbf s\\|^2+\epsilon}\mathbf s } \\]

This is essentially a stochastic, normalized least-squares correction. It costs \\(O(N)\\), and requires no matrix inversion, neural network, or additional energy evaluation.

The interpretation is simple:

- The optimizer predicts the energy change from its current gradient estimate.
- It evaluates one trial configuration.
- The discrepancy corrects the estimated gradient along the tested direction.
- Even an uphill trial provides useful information.

Successful and rejected trials both constrain the local downhill direction. A rejected trial is not wasted information.

There are, of course, qualifications. This formula assumes the energy difference is reasonably approximated by a linear model, and the estimate becomes stale as the configuration moves. I would use forgetting or a limited history, and restrict the update for very large trial displacements.

Furthermore, in a high-dimensional space we cannot reconstruct a full gradient from just a few scalar measurements. We can only learn the gradient components in the subspace we have explored.

But that suggests something more interesting.

### C. Learn the potential only in the two- or three-dimensional subspace we are actually exploring

Suppose the relevant directions are

\\[ \mathbf u_1=\widehat{\mathbf p}, \qquad \mathbf u_2=\text{a new direction orthogonal to }\mathbf u_1. \\]

We can approximately model the energy in this subspace by

\\[ E(a,b)\approx E_0+ g_1a+g_2b+ \frac12\left(H\_{11}a^2+2H\_{12}ab+H\_{22}b^2\right). \\]

That requires learning only five coefficients:

\\[ (g_1,g_2,H\_{11},H\_{12},H\_{22}). \\]

Each attempted move, including failures, supplies one observation. We can update these coefficients incrementally using recursive least squares or a small sliding-window fit.

Now the optimizer can learn not only a downhill direction, but also the orientation of a local narrow valley.

After making progress, it can replace one of its basis vectors and continue learning.

This is a local, moving, low-dimensional surrogate of the potential, rather than a global approximation.

It is attractive because the number of model parameters does not grow with the number of atoms. A 3D quadratic model would require only nine coefficients.

I see two implementation difficulties:

First, we must not combine distant energy samples as if the same quadratic approximation remained valid everywhere. The trust region and forgetting scheme are essential.

Second, an estimated quadratic minimum may lie far outside the explored region or may not exist if the fitted Hessian is indefinite. The predicted move should therefore be limited by a trust radius, with actual energy improvement used to validate it.

I would start with the linear gradient estimator above before building the quadratic model, because the simpler method already tests whether failed moves provide substantial predictive information.

## 5. How I would train such an optimizer

I would initially avoid training the optimizer using backpropagation through the complete optimization trajectory.

Instead, define an algorithm with perhaps 5–20 trainable scalar parameters and optimize those parameters on many example problems.

For example, parameterize the momentum coefficient, exploration amplitude, step-size adaptation and curvature response by a small vector \\(\theta\\).

Then define the meta-objective

\\[ J(\theta)= \mathbb E\_{V,\mathbf x_0} \left[ C\_{\rm relax}(V,\mathbf x_0;\theta) \right]. \\]

Here \\(C\_{\rm relax}\\) could be the number of force evaluations needed to reach a prescribed convergence threshold, with a penalty for failures.

This can be optimized offline using evolutionary strategies, CMA-ES, or other black-box methods. We do not need derivatives of the optimization algorithm, or derivatives of the energy beyond what the inner optimizer normally uses.

And if the parametric rule is not expressive enough, we could subsequently replace the scalar mapping with a tiny MLP.

### A particularly important distinction about training data

You suggested training using terabytes of MD trajectories and forces. Those would certainly be useful, but I would not treat ordinary MD trajectories as the primary training examples.

MD explores phase space according to dynamics; minimization trajectories follow a very different distribution. What constitutes a good optimization step cannot generally be inferred from the force at one configuration.

For example, an MD trajectory might sample a narrow valley and show the local restoring forces, but it doesn't directly tell us whether a proposed large optimization step would have accelerated convergence or sent us into another basin.

I would instead generate training data through optimization rollouts with the actual potential.

The data would include initial configurations, all attempted steps, energies, forces where available, and final convergence outcomes. MD configurations would provide a rich distribution of starting states.

Given access to fast batched force-field evaluations, this becomes a natural large-scale offline training problem.

### Train on a hierarchy of surfaces

[Bohrium | AI for Science with Global Scientists](https://images.openai.com/static-rsc-4/y1LOf_OiMLrH3EfPy0Ea7V-DyDReVdq9svMTNMisGdmEn5XvizElQdWsdDmil1Fynnvrv_NSuYyOinPs9C-ssJdMo6kcVsYNVJqwvsWtZPhaqetULrWyUspP1CZp-sfGo4chRg-v0OUmXH7WNkex8xo8tnIyv9NoXk3AMiAukJM?purpose=inline)

Stage 1: Synthetic landscapes

Rosenbrock valleys, rotated anisotropic quadratics, coupled nonlinear coordinates, periodic potentials, curved manifolds and multimodal functions. Randomize curvature, stiffness ratios, noise, length scales and dimensionality.

These allow millions of cheap training episodes and expose specific weaknesses.

[Tour de Code | Visualization of molecules at surfaces](https://images.openai.com/static-rsc-4/76pt1roChWZGYpBNN6YxygyIy54WK_FDjxuKcaaR6agvyyp_hUj5hjb19AXyKITlUyQUxROQUo_FQHPPVPhQDt9rz9Hhe1Ax6LZrQtcO3zqd0yw2Q-wlNMur11eQeU4a9FlZ2C458nuTk41S8TS46r_laxAPxXkAvPp3VIMGea8?purpose=inline)

Stage 2: Physical potential-energy surfaces

Atomistic force fields with different bonded interactions, intermolecular potentials, adsorbates, flexible molecules and molecular assemblies.

These teach the optimizer the distributions of stiffness, coupling, curvature and force magnitudes found in the intended applications.

[Atomic force microscope reveals chemical ghosts - Berkeley News](https://images.openai.com/static-rsc-4/hCxUAgeYY6jetqv2E2KA2bADTdZua5z5ShiTw3XdAvvLushFqByQlBVkh-jEkNak2Jnp1vzUba_jO5arNIx6kgXAq7gp5BPqRoKcubeHabpLWawWb8r0bpJeedD33YJaS2KPM3Kpr463FWirulWwdkgTDq-Y6AurY1YW0nrWoCA?purpose=inline)

Stage 3: Application-specific optimization

Relaxation of organic molecules on surfaces, molecular assembly and rigid-body manipulation. Sample initial states from realistic MD and Monte Carlo trajectories.

Here I would specialize the previously learned optimizer to the actual force-field family while preserving its generally useful behavior.

One important precaution: the final test set should contain different molecular systems and potential parameters, not merely different initial configurations of the same molecules.

Otherwise, we risk learning an optimization algorithm that works on a small collection of training landscapes but has little broader value.

## 6. Pretraining followed by rapid online adaptation

Your ideal combination of a pretrained prior and rapid adaptation to a particular potential has another direct precedent in the literature: [Learned Optimizers That Continually Adapt](https://opt-ml.org/papers/2023/paper109.pdf?utm_source=chatgpt.com) (OPT 2023 workshop).

That work explores continual adaptation of a pretrained optimizer during its deployment, using a safeguard that falls back on a conventional optimizer when performance deteriorates. It is conceptually relevant, although the adaptation procedure is more expensive than I would want in a lightweight molecular minimizer.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://opt-ml.org\&sz=32)

opt-ml.org



For our case, I would propose something more economical.

The offline training determines the general strategy: how strongly to trust momentum, when to explore, how to respond to force reversal, and how aggressively to increase the step size.

The online learning estimates properties of the particular potential: local curvature, productive directions and characteristic length scales.

This can be implemented with a few state variables, or an online least-squares model, without changing any neural-network weights.

An additional possibility is to give the optimizer a few fast adaptive parameters (perhaps three to five) and let them change according to learned rules. But I would only introduce this if simpler state adaptation proves insufficient.

## 7. What could go wrong, and how to test it?

There are several limitations worth considering before investing heavily in training.

A learned optimizer need not beat a good L-BFGS implementation. In one systematic atomistic benchmark, L-BFGS outperformed CG and FIRE for local optimization, often by substantial margins. FIRE 2.0 has demonstrated advantages over CG on other problems, so performance is strongly problem-dependent.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://pubs.acs.org\&sz=32)

ACS Publications+1



The useful target is not to defeat BFGS in its ideal quadratic regime. It is to improve efficiency during difficult, strongly nonlinear parts of a trajectory without sacrificing reliable final convergence.

A large learned optimizer may not pay for itself. A particularly relevant 2025 study by Trukhan, Mazhnik and Oganov uses deep reinforcement learning for crystal relaxation. The method can compete with classical optimizers in some regimes, but generalization and strict convergence near minima remain challenging. The authors also suggest switching to conventional optimizers near convergence.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.nature.com\&sz=32)

Nature



Their approach employs graph neural networks, which is substantially more complicated than the history-based controller proposed here.

Global optimization is a different problem. If we always reject uphill steps, our optimizer cannot escape a strict local minimum, regardless of how much it learns. To explore different basins, we would need some mechanism permitting uphill transitions, restarts or sufficiently nonlocal proposals. I would keep this separate from the initial problem of navigating curved valleys efficiently.

The computational environment matters. On a GPU, even a tiny MLP evaluated separately for every atom could create unwanted memory traffic and kernel-launch overhead. A single scalar controller evaluated once per configuration, followed by a fused vector-update kernel, would be much more attractive.

For initial benchmarks, I would use the following comparison:

| Optimizer                               | Energy only | Forces | Online state                              |
| --------------------------------------- | ----------- | ------ | ----------------------------------------- |
| Your momentum + random-search heuristic | Yes         | No     | Successful-step history                   |
| Learned momentum + gradient estimator   | Yes         | No     | Momentum and estimated gradient           |
| Low-dimensional quadratic surrogate     | Yes         | No     | Local model coefficients                  |
| (1+1)-CMA-ES                            | Yes         | No     | Evolution path, covariance                |
| FIRE 2.0                                | Optional    | Yes    | Velocity, timestep                        |
| L-BFGS                                  | Optional    | Yes    | Limited secant history                    |
| Learned FIRE controller                 | Optional    | Yes    | Momentum, force history, learned controls |

The crucial performance metrics would be total potential evaluations, actual wall time, convergence success rate and final energy, preferably reported as distributions across many independent problems. For force-based optimization, the final maximum force should also be measured.

In particular, for your GridFF-type applications where analytic forces are inexpensive, the force-based version may have a substantial advantage. A single force evaluation gives an \\(N\\)-dimensional vector of information, whereas an energy-only evaluation provides one scalar. The latter becomes most attractive when forces are unavailable, substantially more expensive, or when the true search space is only a few-dimensional set of collective coordinates.

## 8. Which direction I would investigate first

I see three promising candidate algorithms:

Learned FIRE

Keep the FIRE philosophy, but learn momentum and step-size control from force-history features. Very cheap inference, straightforward to compare against FIRE 2.0.

Gradient-free momentum with online gradient estimation

The most direct extension of your previous heuristic. Exploits information from rejected moves without building a full Hessian or surrogate. Particularly elegant in low-dimensional search spaces.

Learned low-rank local model

Track a few useful collective directions and approximate energy or force behavior in their subspace. More sophisticated, but potentially very effective in the curved-canyon problem.

Of these, I find the third direction scientifically the most interesting. It is a useful middle ground between inertial motion, local surrogate modeling, and quasi-Newton methods, without committing to a large Hessian approximation.

But I would implement the first two before training any neural network. Parameterize them with a handful of coefficients, optimize those coefficients offline over thousands of synthetic and physical landscapes, and see how far such a simple approach can go.

Only after reaching the limits of fixed coefficients would I introduce a tiny MLP or symbolic program search. At that point we would know exactly what additional information and decision-making capability the learned optimizer needs.

The research opportunity I find most compelling is to discover whether an optimizer with just a few history vectors and perhaps a dozen learned scalars can consistently outperform its handcrafted equivalents on a well-defined family of physical potential-energy surfaces. That would be both useful scientific software and a more interpretable result than a large neural network trained to output atomic displacements.

---
---

# USER


OK, good research

I kida like most the frameworks of Lion
[https://arxiv.org/pdf/2302.06675](https://arxiv.org/pdf/2302.06675)
baisucally find alternative formulations somilar to exsiting stimulated anealing, Adam, or FIRE which are simple mathematical expresion with few learned paramters which have some physical meaning (like interia), rhat then using opaque and messy MLP to guess next move. It still make sense to me to use MLP-controler to guess how to adjust the paremters as we move long the trajectroy to react to situation (are we stuck, are we in long narrow valye, are we near minimum, how much is the valy curved perpendicular to softest direction ... ) &#x20;

ad "There is a trade-off: a few global scalars cannot adapt separate regions of a molecule with very different stiffness. One possible extension is to run the same controller independently for different blocks: rigid-body translation, rotation, internal molecular coordinates, substrate degrees of freedom, etc.:

I would stress the basic philosophy - we do not want to tran anything with would try tu undetand the underlaying problem or memorize the potentil. The controler should not try to undestand molecular geometry/topology (it would make it too complicated, dostly and no transferable - this would be as costly as traning machine learning potential), we want minimalistic controler which is just a litbnle bit smarter and more complex then FIRE,ADAM,POWER or Lion, it shold not see the underlaying problem (molecular geometry and topology) at all, it should decide purely based on dynamical variables at hand (e.g. some memory of dot products, parallel and peprpedictualt projections of some vectors - displacement, velocity, curvature along the step etc. )

---
---

# Chat GPT 6

I agree with this philosophy, and I think it gives us a much more precise definition of the project:

Discover simple, physically interpretable optimization algorithms that adapt to the local dynamics of an arbitrary energy landscape, without learning anything about the underlying system.

The optimizer should be completely agnostic to whether it is relaxing a molecule, minimizing an electronic structure, fitting parameters, or optimizing a trajectory. It sees only a sequence of positions, energies, forces (when available), attempted displacements and its own internal state.

I would distinguish two things:

1. The optimizer itself: a compact mathematical update rule, analogous to FIRE, ADAM, Lion or simulated annealing, with perhaps 5–15 meaningful parameters.
2. The adaptation mechanism: a small controller that observes the history of optimization dynamics and adjusts those parameters. Initially, this could be handcrafted or symbolically discovered; a tiny MLP is an optional later extension.

There is an important lesson from Lion: its main attraction for us is not necessarily the `sign()` operation, but the discovery methodology. Symbolic program search found a useful algorithm and simplified it into a few mathematical operations.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://papers.neurips.cc\&sz=32)

NeurIPS Papers+1



In fact, I would explicitly avoid Lion's coordinate-wise `sign()` as a default operation for physical minimization. It is not rotationally equivariant: rotating Cartesian coordinates changes the physical direction of the update.

Instead, our building blocks should be vector operations with clear geometrical meanings: scalar products, norms, projections, normalization, linear combinations, momentum, and changes in force.

This suggests a particularly compact architecture.

## A concrete minimalist design

### 1. The optimizer should operate in a geometrical basis defined by its own motion

Suppose we maintain only:

- Current force \\(\mathbf F_t=-\nabla E(\mathbf x_t)\\).
- Previous forces and displacements.
- Momentum \\(\mathbf m_t\\), representing the history of useful movement.
- A characteristic step length \\(h_t\\).
- A few exponentially averaged scalar diagnostics.

The current direction of travel is

\\[ \mathbf u_t=\frac{\mathbf m_t}{\\|\mathbf m_t\\|}. \\]

We decompose the current force into

\\[ \mathbf F\_\parallel=(\mathbf F\cdot\mathbf u)\mathbf u, \qquad \mathbf F\_\perp=\mathbf F-\mathbf F\_\parallel. \\]

Now consider this candidate update:

\\[ \boxed{ \Delta\mathbf x_t = a_t\mathbf m_t+ \frac{h_t}{\bar F_t} \left( b_t\mathbf F\_\parallel+ c_t\mathbf F\_\perp \right) } \\]

where \\(\bar F_t\\) is a running force scale and \\(a_t,b_t,c_t\\) are dimensionless coefficients.

This is already more expressive than simple momentum mixing, while retaining a transparent physical interpretation.

| Parameter | Physical meaning                                                                |
| --------- | ------------------------------------------------------------------------------- |
| \\(a_t\\) | Inertia: how much previously established movement to preserve                   |
| \\(b_t\\) | Response to the force along the direction of travel; acceleration or braking    |
| \\(c_t\\) | Response to the transverse force; steering and correction of lateral deviations |
| \\(h_t\\) | Characteristic exploration length                                               |

If \\(\mathbf m=0\\), we simply initialize the direction from the force.

A single parameter \\(c_t\\) is particularly interesting. In a narrow valley, an instantaneous transverse force may indicate either a restoring force from a stiff valley wall or a necessary correction to follow a bending valley.

Those two situations require different responses.

### 2. The most useful signals may come from temporal correlations, not instantaneous forces

This is where I think a slightly smarter FIRE could offer real benefits.

Imagine two scenarios.

Scenario A: oscillating between the walls of a straight, narrow valley.

The transverse force tends to reverse repeatedly. The appropriate response is to damp the unwanted oscillation rather than follow every reversal.

Scenario B: following a curved valley.

The transverse correction can remain consistently directed toward the turning valley over several steps. We may want to retain momentum while continuously steering.

This suggests measuring the persistence of the transverse force:

\\[ C\_\perp= \frac{ \mathbf F\_{\perp,t}\cdot\mathbf F\_{\perp,t-1} }{ \\|\mathbf F\_{\perp,t}\\| \\|\mathbf F\_{\perp,t-1}\\|+\epsilon }. \\]

Both forces should be projected onto the transverse subspace defined by the current direction \\(\mathbf u_t\\), to make the comparison consistent.

Then:

- \\(C\_\perp\approx-1\\): alternating transverse forces, suggestive of oscillation.
- \\(C\_\perp\approx+1\\): persistent transverse force, potentially indicating the need to steer.
- \\(C\_\perp\approx0\\): little directional persistence.

Of course, this distinction isn't definitive. A persistent transverse force could also arise from motion in a consistently wrong direction. The controller would have to combine this signal with energy progress and force alignment.

But importantly, we are no longer asking a neural network to understand valley geometry. We are asking it to learn how to respond to a measurable temporal pattern.

### 3. A small collection of useful dynamical observables

I would start with approximately six scalar features.

| Observable                   | Formula or definition                                                | What it tells us                                          |
| ---------------------------- | -------------------------------------------------------------------- | --------------------------------------------------------- |
| Force–motion alignment       | \\(\cos(\mathbf F,\mathbf m)\\)                                      | Moving downhill, sideways or uphill?                      |
| Force-direction persistence  | \\(\cos(\mathbf F_t,\mathbf F\_{t-1})\\)                             | Is the driving direction stable or rapidly changing?      |
| Transverse-force persistence | \\(C\_\perp\\)                                                       | Steering versus oscillation                               |
| Longitudinal curvature       | \\(-\Delta\mathbf F\cdot\Delta\mathbf x/\\\|\Delta\mathbf x\\\|^2\\) | Local stiffness along the last displacement               |
| Relative force magnitude     | \\(\\\|\mathbf F_t\\\|/\overline{\\\|\mathbf F\\\|}\\)               | Is the force becoming unusually weak or strong?           |
| Energy-prediction accuracy   | Actual versus force-predicted \\(\Delta E\\)                         | Is the previous step exceeding the locally linear regime? |

Curvature must be normalized by suitable running force and displacement scales before entering a transferable controller.

There is also useful information in the transverse component of \\(\Delta\mathbf F\\). Under the local approximation

\\[ \Delta\mathbf F\simeq-H\Delta\mathbf x, \\]

it measures coupling between the direction of displacement and other directions. However, it does not directly measure the transverse stiffness or the curvature of the valley floor.

That distinction is important. With only a history of displacements in one direction, we cannot reconstruct the Hessian in directions we have never explored. To learn transverse stiffness, the optimizer must have experienced some transverse displacement.

Similarly, a small force alone cannot distinguish proximity to a minimum from a flat region far away. Some situations are intrinsically ambiguous given the available observations.

These are limitations of the information, not of the controller architecture.

## 4. What the tiny MLP would actually do

The MLP would never output \\(\Delta\mathbf x\\) or any vector.

It would receive a handful of normalized scalar diagnostics and output adjustments to the coefficients of our mathematical update.

For instance,

\\[ (a_t,b_t,c_t,\delta h_t) =\pi\_\theta(\mathbf z_t) \\]

where \\(\mathbf z_t\\) contains the six dynamical observables.

The network could have just one hidden layer of eight neurons. That is roughly 90 parameters, evaluated once per step, independently of the number of optimization coordinates.

But I would go further: rather than allowing the MLP to choose arbitrary coefficients, let it make small corrections to reasonable baseline values:

\\[ a_t=a_0+\delta a_t,\qquad b_t=b_0+\delta b_t,\qquad c_t=c_0+\delta c_t. \\]

We can bound the corrections to avoid pathological behavior.

This gives us a graceful progression from a fixed heuristic to an adaptive one.

An even better possibility is to use symbolic regression on the trained policy to produce explicit equations such as

\\[ c_t= c_0+ c_1\max(0,C\_\perp) + c_2\max(0,-C\_\perp). \\]

This is only an illustrative expression; the actual coefficients and functional form would be discovered.

The interesting result would be if such an explicit controller consistently outperformed a fixed-parameter FIRE-like algorithm.

There is precedent for this optimism: [Reverse engineering learned optimizers reveals known and novel mechanisms](https://proceedings.neurips.cc/paper_files/paper/2021/hash/a57ecd54d4df7d999bd9c5e3b973ec75-Abstract.html?utm_source=chatgpt.com) (NeurIPS 2021) found that learned optimizers often internally develop recognizable mechanisms such as momentum, gradient clipping, and adaptive learning-rate schedules. We could build those interpretable mechanisms into the search space from the beginning.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://proceedings.neurips.cc\&sz=32)

NeurIPS Proceedings



## 5. I would actually explore two complementary discovery approaches

Approach A — Parametric physics-inspired optimizer

Start with the explicit equation above. Train its coefficients and possibly a tiny controller to optimize performance across different landscapes.

This has a strongly constrained search space and is easy to debug.

Approach B — Symbolic algorithm discovery, as in Lion

Instead of assuming the form of the update, define a small vocabulary:

\\[ \begin{gathered} +,\\;-,\\;\times,\\;/,\\;\min,\\;\max,\\\ \text{dot},\\;\text{norm},\\;\text{normalize},\\\ \text{project},\\;\text{EMA},\\;\text{clip}. \end{gathered} \\]

Search over short programs using these operations, with access only to the dynamical state.

Optimize not only convergence performance but also expression complexity. A shorter program with slightly worse convergence could be preferable to a complicated one.

The symbolic approach is intellectually more ambitious and could discover rules we have not anticipated. But I would constrain it to rotationally equivariant vector expressions and dimensionless scalar diagnostics.

The output should also respect basic physical scaling. For example, multiplying the potential energy by a constant should not fundamentally change the optimization trajectory, provided the optimizer adjusts the characteristic force scale appropriately.

There is a more recent paper on exactly this general research direction: [Symbolic Discovery of Iterative Algorithms: A Continuous Latent Space Bayesian Optimization Framework](https://arxiv.org/abs/2607.01552?utm_source=chatgpt.com) (Liu and Mitrai, July 2026). It searches for iterative update rules in symbolic form using a learned latent representation of programs. Its discovery machinery is relatively elaborate, but the resulting optimizer is still a mathematical expression.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arXiv



## 6. One design decision I think is especially important

I would not start by letting the MLP control too many independent parameters.

Instead, test whether a single extra degree of freedom beyond ordinary FIRE already provides a significant improvement.

For example, keep conventional momentum and timestep adaptation, but allow the ratio

\\[ r_t=\frac{c_t}{b_t} \\]

between transverse and longitudinal force response to vary dynamically.

Does this allow the optimizer to retain its forward momentum in long valleys while suppressing transverse oscillations?

If so, we would have identified a genuinely useful dynamical mechanism, rather than just fitted a complicated optimizer.

Then gradually introduce other degrees of freedom and test their marginal contributions.

This would also help determine whether learning itself is necessary. Perhaps a fixed symbolic formula using a couple of observables is enough.

And the same philosophy extends to your energy-only heuristic: its learned rule should operate on accepted/rejected displacement histories, energy changes, directional correlations and estimated local response, not on any representation of the physical system.

The central research question, in my view, is which small collection of dynamical observables contains enough information to make better optimization decisions than FIRE, ADAM or simple stochastic search. Once those observables are identified, learning their optimal combination may be relatively easy.

That feels like the right level of abstraction: not machine learning of the potential-energy surface, but machine learning of a compact control law for optimization dynamics.


---
---

# USER

OK, I agree with this approach, but the explicit formula you sketched seesm to me too simple constrained and oo similar to exsting algorthms, I do not see much space for improvement.

I should also clarify one thing - the controler should look on trajectroy within reference frame of tis one walk and pose (like current step direction, previous step direction define coordinate system sub-sapace). But for numerical efficienty would shoe avoid storying whole stat-vectors if we optimize e.g. molecule with 1000 atoms or machine-learning-model with 1 million parameters we do not want to store all these state fectos for 100-1000 iterations. We want to store only limited information.

I think we should disitnguies explictly the Gradietn-free and graidient-full version. Gradietn-full version will be perhaps no much different from FIRE (Fire seems to be already close to optimal in many cases, what make sense is somehow selectively re-weight different DOFs to compensate non-equal curvatures, but that require some insight into structure of particualr problem - repsp. al least tracking history of these dofs independnetly)

For gradient-less optimization there is much more room for heuristics I thin as it seems to me as unsolved problem and most gradient-less methods are still very bad.

We should aim for graeidnt-free method which aproach efficient of gradient full method - that would provide huge advantage as gradients are often difficutl to impleemnt analytically and in any case require aditional computing cost. 

I was alsways dereaming aporximating gradients using previoous samples, which more or less the inertia does. It is ease to know gradient component (projection) along the last step, simply <F|dij>=(Ei-Ej)/|ri - rj|, but that is not the full story. It make sesne to try to exoplere around to estimate gradeint and curvature in different direction. It would be good it the model bild some small cheap low-dimentional representation of gradient and curvature (hessian). Perhaps we can consider few (2-3) soft modes (vally direction and in which direction the valye is recurved) and the rest can be seen as some "hard elastic bath" (the other stiff modes. The stiffnes of soft modes we estimate from our path, and the stiffnes of hard modes (then steep slops aroundnarrow valye) we estimate from errorous samples in radom direction roughly perpendicult to the path. when we hit a wall (the valye turn sharply) we ne to spend sampled on searching new soft descent direction. This gives tradeoff between exploration and exploitation.

 Natureally when descent is eficient we are in explitation phase we want to go down the hill as fast as possible. When we get stuck and energy stop decreasing and start increasing we must search around in different direction. (it is like heruristic of stimulated anealing by in spce of descent direction rather than positions).

---
---

# Chat GPT 6

I think this is a substantially more interesting direction than the adaptive-FIRE controller we discussed earlier.

The central idea I take from your proposal is:

A gradient-free optimizer should learn a moving, low-dimensional representation of the local descent landscape, using its own trajectory and the information gained from both successful and unsuccessful moves.

Instead of approximating the entire gradient or Hessian, it should identify a few directions in which useful motion occurs, estimate how those directions evolve, and treat the remaining degrees of freedom statistically.

This naturally produces the exploration–exploitation trade-off you describe: exploit a known downhill corridor while it remains predictive; spend additional energy evaluations identifying new directions when the model fails.

There is a particularly important distinction here: we should learn both the direction of travel and how the useful subspace itself moves through configuration space. That is more expressive than ordinary momentum, and potentially much more efficient than random search.

## 1. First, a limitation worth making explicit

Your objective of approaching gradient-based efficiency using energy evaluations alone is plausible for some classes of landscapes, but impossible uniformly for arbitrary high-dimensional functions.

Consider a locally linear potential,

\\[ E(\mathbf x)=E_0+\mathbf g\cdot\mathbf x, \\]

where \\(\mathbf g\\) is an arbitrary unknown \\(N\\)-dimensional gradient.

Each energy difference gives one scalar measurement,

\\[ \Delta E=\mathbf g\cdot\Delta\mathbf x. \\]

With \\(k\\) independent measurements we learn at most \\(k\\) projections of \\(\mathbf g\\). The other \\(N-k\\) components remain unknown.

For an isotropically distributed unknown gradient, \\(k\\) random orthogonal measurements capture, on average,

\\[ \frac{\\|\mathbf g\_{\rm estimated}\\|^2}{\\|\mathbf g\\|^2} =\frac{k}{N}. \\]

For a million-dimensional optimization problem, three random measurements reveal almost nothing about an arbitrary gradient.

This is a fundamental information limitation, not a weakness that a learned heuristic can eliminate. The dimensional dependence of zeroth-order optimization is well established, although it can be significantly reduced for structured problems.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://ieeexplore.ieee.org\&sz=32)

IEEE Xplore+1



However, your curved-valley example provides exactly the structure we could exploit.

If the optimizer already knows a useful direction and the next useful direction is strongly correlated with it, each new evaluation only needs to correct its existing prediction, rather than reconstruct an entirely new gradient.

That is where trajectory memory can provide a major advantage.

## 2. Your soft-mode plus elastic-bath model is a good starting point

Let the optimization problem have \\(N\\) coordinates, but suppose that locally the interesting motion occurs in a subspace of dimension \\(k=2\\) or \\(3\\).

Maintain an orthonormal basis

\\[ U=[\mathbf u_1,\ldots,\mathbf u_k], \qquad U^\mathsf TU=I. \\]

The directions could have the following dynamical interpretation:

- \\(\mathbf u_1\\): recent successful direction of descent.
- \\(\mathbf u_2\\): direction in which the trajectory bends, approximately the normal to the valley.
- \\(\mathbf u_3\\): a newly explored direction, possibly replacing an older one.

All other degrees of freedom belong to the complementary subspace.

We could approximate the local energy by

\\[ \boxed{ E(\mathbf x+\mathbf s)\approx E_0+ \mathbf g_s^\mathsf T\mathbf a+ \frac12\mathbf a^\mathsf TH_s\mathbf a+ \frac12 k_b\\|\mathbf r\\|^2 } \\]

where

\\[ \mathbf a=U^\mathsf T\mathbf s, \qquad \mathbf r=(I-UU^\mathsf T)\mathbf s. \\]

Here:

- \\(\mathbf g_s\\) is the estimated gradient within the soft subspace.
- \\(H_s\\) is a \\(2\times2\\) or \\(3\times3\\) curvature matrix.
- \\(k_b\\) is the effective curvature of the hard bath.

The model assumes the transverse gradient is negligible near the valley floor. We should relax that assumption when exploring unfamiliar directions.

This is not a representation of the molecular structure or the actual Hessian. It is a compact model of how energy responds to motion in the directions that the optimizer has discovered.

For a two-dimensional soft subspace, it contains only two gradient coefficients, three Hessian coefficients, and one bath stiffness.

Six scalar parameters describing the local landscape.

The basis vectors are the expensive part, but we only store two or three, not hundreds of historical configurations.

Approximate memory for three float32 direction vectors

| 3,000 coordinates (1,000 atoms) | 36 KB |
| ------------------------------- | ----- |
| 1,000,000 coordinates           | 12 MB |

Excludes current coordinates and temporary workspace. Storage is \\(O(kN)\\), independent of trajectory length.

This is similar in spirit to limited-memory optimizers, but our basis is specifically intended to represent a moving valley and its local curvature, not just accumulate descent history.

## 3. I would go beyond a quadratic model: explicitly represent a curved valley

This is perhaps the most interesting extension of your idea.

A quadratic approximation is a natural starting point, but it cannot describe a curved canyon over a significant distance. Its valley is always straight.

Yet we can represent a curved valley using only a few more parameters.

Consider coordinates \\(s,z\\) defined relative to the current direction of travel and its estimated normal:

\\[ \mathbf x(s,z)=\mathbf x_0+s\mathbf t+z\mathbf n. \\]

Imagine the minimum-energy valley follows approximately

\\[ z\_{\rm valley}(s)=\frac{K}{2}s^2 \\]

where \\(K\\) describes its geometrical curvature.

Now model the energy as

\\[ \boxed{ E(s,z)=E_0+ g_s s+\frac12 k_s s^2+ \frac12 k\_\perp \left(z-\frac K2 s^2\right)^2 } \\]

This expression is simple, but it already describes a rather nonlinear landscape.

Schematic curved valley in a moving local coordinate system. The energy increases steeply when we deviate from the curved centerline.

The four parameters have straightforward meanings:

| Parameter      | Meaning                                    |
| -------------- | ------------------------------------------ |
| \\(g_s\\)      | Downhill slope along the valley            |
| \\(k_s\\)      | Change of slope along the valley           |
| \\(k\_\perp\\) | Transverse stiffness of the canyon walls   |
| \\(K\\)        | Geometrical curvature of the valley's path |

Notice the distinction between energy curvature \\(k_s,k\_\perp\\) and geometrical curvature \\(K\\). These are physically different quantities.

For example, a valley can be very flat along its centerline (\\(k_s\\) small), extremely stiff transversely (\\(k\_\perp\\) large), and strongly curved (\\(K\\) large).

This is exactly the difficult situation you describe.

### Why predicting the bend matters so much

Suppose we take a straight step of length \\(\ell\\) along the current tangent.

The departure from the valley centerline is approximately

\\[ \delta z\approx\frac{K\ell^2}{2}. \\]

The unwanted transverse energy penalty is therefore

\\[ \boxed{ \Delta E\_{\rm wall} \approx\frac{k\_\perp K^2\ell^4}{8}. } \\]

It grows as the fourth power of step length.

So a modest increase in step length can be disastrous, even when the downhill slope is nearly unchanged.

But if we predict the valley curvature and choose

\\[ \Delta\mathbf x =\ell\mathbf t+\frac{K\ell^2}{2}\mathbf n, \\]

we can eliminate this leading transverse error in the idealized model.

We do not need the full Hessian to do that. We only need a reasonable estimate of the turning direction and its curvature.

In a smooth valley, the remaining geometric prediction error is then higher order in \\(\ell\\). Whether this actually improves efficiency depends on how reliably we estimate \\(K\\) from imperfect trajectories.

### How do we estimate the bending without gradients?

From the geometry of previously successful displacements.

Define the normalized directions of two consecutive accepted steps:

\\[ \mathbf t\_{t-1} =\frac{\Delta\mathbf x\_{t-1}}{\\|\Delta\mathbf x\_{t-1}\\|}, \qquad \mathbf t_t =\frac{\Delta\mathbf x_t}{\\|\Delta\mathbf x_t\\|}. \\]

Their turning angle is

\\[ \theta_t=\arccos(\mathbf t_t\cdot\mathbf t\_{t-1}). \\]

For approximately uniform steps along a smooth path,

\\[ K_t\approx \frac{\theta_t} {(\ell_t+\ell\_{t-1})/2}. \\]

The turning normal is obtained by projecting the change in tangent perpendicular to the current tangent.

We have therefore estimated the direction and curvature of the path from just two consecutive displacement vectors.

I would smooth this estimate over several steps using a recurrent state, because zigzagging across a canyon can otherwise be mistaken for genuine valley curvature.

The most important caveat: a trajectory of accepted downhill steps is not necessarily the actual valley centerline. The curvature estimate can be contaminated by transverse relaxation. That is where the history of rejected steps, and occasional transverse probes, becomes valuable.

## 4. Failed moves should be used as measurements of the canyon

We now have a local model that predicts both a downhill direction and its likely bending.

What happens when a predicted move increases energy?

We should try to distinguish three possibilities:

| Observation                                    | Possible explanation                                | Useful response                                        |
| ---------------------------------------------- | --------------------------------------------------- | ------------------------------------------------------ |
| Energy rises sharply even for a small step     | Stiff transverse motion, wrong direction            | Improve steering; estimate wall curvature              |
| Energy rises only for longer steps             | Overshoot, changing curvature, nonlinear energy     | Shorten step or improve path extrapolation             |
| Energy decreases, but much less than predicted | Model becoming stale, approaching a turn or minimum | Update local model; possibly explore another direction |

A failed move should not automatically destroy the useful momentum.

Instead, it supplies an observation of the form

\\[ \Delta E_i =E(\mathbf x+\mathbf s_i)-E(\mathbf x). \\]

We can fit our four-parameter curved-valley model, or the more general low-dimensional quadratic model, to these observations.

Crucially, we do not need to store all historical displacements.

For a fixed local coordinate frame, the quadratic model is linear in its coefficients. Each new observation contributes a small regression row such as

\\[ \boldsymbol\phi_i= \left( s_i,\ z_i,\ \frac12s_i^2,\ s_i z_i,\ \frac12z_i^2,\ \frac12\\|\mathbf r_i\\|^2 \right). \\]

We can maintain exponentially discounted normal-equation statistics,

\\[ A\leftarrow\lambda A+\boldsymbol\phi_i\boldsymbol\phi_i^\mathsf T, \qquad \mathbf b\leftarrow\lambda\mathbf b+\boldsymbol\phi_i\Delta E_i. \\]

Solving \\(A\boldsymbol\theta=\mathbf b\\), with regularization, recovers the local model coefficients. Numerically, a small QR or Cholesky implementation with conditioning checks would be preferable to explicitly inverting \\(A\\).

For six model coefficients, \\(A\\) is just \\(6\times6\\).

This gives an inexpensive way to learn from dozens of evaluations while storing only a handful of scalars.

There is an important complication: these statistics describe one local coordinate frame and one reference point. As the optimizer moves and rotates its basis, observations cannot simply be retained unchanged. We can translate the model within its validity region, but when the frame changes substantially we should discount or reset the fit. Alternatively, retain a small cache of recent full displacement vectors. The memory would still be \\(O(kN)\\), not \\(O(TN)\\).

For the curved-valley model, the fitting problem is nonlinear in \\(K\\) and \\(k\_\perp\\). A practical solution is to estimate the geometric curvature \\(K\\) from trajectory history, then fit the remaining energy parameters with a small regression.

### Measuring directional gradient and curvature separately

There is another useful trick that doesn't require a surrogate model.

At a fixed position, evaluate two opposite probes:

\\[ E\_+=E(\mathbf x+\delta\mathbf u),\qquad E\_-=E(\mathbf x-\delta\mathbf u). \\]

Then

\\[ \boxed{ g_u\approx\frac{E\_+-E\_-}{2\delta}, \qquad k_u\approx\frac{E\_++E\_--2E_0}{\delta^2}. } \\]

Two evaluations give the directional gradient and curvature, up to finite-difference error.

That costs more than one random trial, but gives substantially cleaner information. In particular, it separates an uphill slope from transverse stiffness.

One could learn when an additional opposite-direction probe is worth its cost.

This is a natural exploration decision for a controller.

## 5. Exploration should search for new directions, not necessarily new positions

I especially like your suggestion of simulated annealing in the space of descent directions.

There are actually two different optimization problems:

\\[ \underbrace{\min\_{\mathbf x}E(\mathbf x)}\_{\text{physical optimization}} \qquad\text{and}\qquad \underbrace{\min\_{\\|\mathbf d\\|=1} E(\mathbf x+\ell\mathbf d)}\_{\text{direction search}}. \\]

During ordinary exploitation, we optimize the physical position using the current best direction.

When progress deteriorates, we can temporarily hold the current best position fixed and explore candidate directions. The direction-search algorithm can accept a worse directional proposal without moving the physical system uphill.

That distinction is useful because we do not have to sacrifice our best configuration to explore different ways forward.

### A concrete proposal mechanism

Let \\(\mathbf t\\) be the current valley tangent, \\(\mathbf n\\) its estimated turning normal, and \\(\mathbf v\\) a random unit vector orthogonal to the currently known soft subspace.

Generate proposals of the form

\\[ \boxed{ \Delta\mathbf x= \ell\mathbf t+ \frac12K\ell^2\mathbf n+ \alpha\mathbf v } \\]

with adaptive \\(\ell\\) and \\(\alpha\\).

The three terms have different meanings: exploitation along the valley, deterministic steering based on its curvature, and exploration of an unknown direction.

We could generalize this to several remembered soft directions.

How the optimizer responds to an energy evaluation is the essential part:

Productive descent

When the energy decreases as predicted, preserve the existing subspace, increase the exploration length cautiously, and extrapolate the valley trajectory. Most iterations should require just one new energy evaluation.

Model correction

When the observed energy differs substantially from the prediction, use the residual to update the local model. Try corrections within the known two- or three-dimensional subspace before introducing a new full-dimensional direction.

Search for a new soft direction

After repeated failures, introduce a random direction \\(\mathbf v\\) orthogonal to the known subspace. Compare energies for positive and negative perturbations when useful. If a new direction proves productive, incorporate it into the soft basis and forget an older, less useful direction.

The resulting algorithm continually alternates between exploiting its existing knowledge and modifying that knowledge.

One complication is that exploration in the unknown complement becomes harder in high dimensions.

If there is only one undiscovered soft direction in a million-dimensional hard bath, random probes will almost never align strongly with it. Their measured curvature will mostly reflect the average stiffness of the other directions:

\\[ \mathbb E\_{\mathbf v} [\mathbf v^\mathsf T H\_\perp\mathbf v] = \frac{\operatorname{tr}H\_\perp}{N-k}, \\]

for uniformly distributed unit vectors in the complementary subspace.

Consequently, the hard-bath approximation is useful for estimating the typical penalty associated with random exploration, but it does not by itself identify the orientation of an exceptional soft mode.

The algorithm will work best when new soft directions emerge smoothly from the existing trajectory, as in a bending valley. Discovering an entirely unrelated soft direction may require many evaluations.

## 6. What should the learned controller decide?

I would now give the controller a different role from what we discussed earlier.

It should not merely adjust momentum or dampening. It should decide how to allocate the next energy evaluation.

For example:

| Controller output      | Meaning                                                               |
| ---------------------- | --------------------------------------------------------------------- |
| \\(\ell\\)             | How far to extrapolate along the current valley                       |
| \\(\alpha\\)           | How much to explore perpendicular to known directions                 |
| \\(\tau\\)             | How strongly to trust and retain the local model                      |
| \\(p\_{\rm explore}\\) | Probability of spending the next evaluation on a new direction        |
| \\(p\_{\rm mirror}\\)  | Whether to evaluate an opposite probe to separate slope and curvature |

The controller's inputs could be progress per evaluation, rejection frequency, angular changes of successive steps, local curvature estimates, model-prediction errors, and uncertainty of the fitted coefficients.

We should not confuse a high rejection rate with a universally bad situation. Near a minimum, many rejected steps are inevitable, while in a sharply turning valley repeated rejections may indicate that the exploration direction needs to change.

I would therefore make model reliability a central input.

For instance, define a trust ratio

\\[ \rho_t= \frac{E_t-E\_{\rm trial}} {E_t-E\_{\rm predicted}}, \\]

when the predicted improvement is sufficiently positive.

If \\(\rho\approx1\\), the local model is working. If \\(\rho\\) is small or negative, the model is becoming unreliable.

This is standard trust-region reasoning, but the choice of which new directions to explore when trust deteriorates could be learned.

One particularly attractive formulation would be a controller that estimates the expected benefit of exploitation versus exploration per additional energy evaluation. That directly optimizes what we actually care about.

## 7. Existing work: this is closer to the literature than it first appears

After looking further into derivative-free optimization, I found several particularly relevant precedents.

[Optimization by Moving Ridge Functions (OMoRF)](https://www.tandfonline.com/doi/full/10.1080/0305215X.2021.1886286?utm_source=chatgpt.com) — 2021

Probably the closest match to your moving low-dimensional landscape idea. It builds local surrogate models in reduced subspaces, updating those subspaces as the trust region moves. It explicitly does not require a globally fixed low-dimensional structure.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.tandfonline.com\&sz=32)

Taylor & Francis Online



The important difference in our proposed approach would be emphasizing the dynamics and geometry of the preceding descent trajectory as the source of useful directions, rather than primarily constructing a statistical reduced model.

[ASEBO — Adaptive ES-Active Subspaces](https://papers.neurips.cc/paper_files/paper/2019/hash/88bade49e98db8790df275fcebb37a13-Abstract.html?utm_source=chatgpt.com) — NeurIPS 2019

Dynamically learns useful subspaces for gradient-free optimization and adjusts the balance between exploration within those subspaces and exploration outside them. This directly addresses the exploration–exploitation issue you raised.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://papers.neurips.cc\&sz=32)

NeurIPS Papers



[Limited-Memory CMA-ES (LM-CMA)](https://arxiv.org/abs/1404.5520?utm_source=chatgpt.com) — Loshchilov, 2014

Stores only a limited number of important direction vectors instead of a full covariance matrix, bringing memory and sampling arithmetic to \\(O(kN)\\). This is especially relevant to your million-parameter example.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://doi.org\&sz=32)

DOI



[Random-subspace Ellipsoidal Trust-region Optimization (RETRO)](https://www.yiwenchen.me/publications/random-subspace-ellipsoidal-trust-region-methods-for-derivative-free-optimization/?utm_source=chatgpt.com) — October 2026

A very recent method that adapts the geometry of the trust region using information about which subspaces have been explored. Its motivation includes an important issue: successful agreement between predicted and actual energy reductions does not necessarily imply that the chosen subspace is well aligned with the true gradient.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.yiwenchen.me\&sz=32)

Yiwen Chen



The basic ingredients of our proposal are therefore not unprecedented. I would not claim novelty before comparing their actual algorithms.

Nevertheless, these methods often use machinery that is more general than what you want. Your approach could be intentionally narrower: a very small stateful optimizer focused on following smooth but highly curved valleys, with perhaps two or three moving directions and a tiny explicit energy model.

## 8. What I would build as the first experiment

I would temporarily set aside the MLP and symbolic regression. First establish whether the low-dimensional moving-valley model has useful predictive power.

Call the experimental method Moving Valley Search for now.

The first implementation could maintain only:

| Current best configuration                     | \\(N\\) coordinates |
| ---------------------------------------------- | ------------------- |
| Valley tangent and normal                      | \\(2N\\) values     |
| Optional new exploration direction             | \\(N\\) values      |
| Local energy/curvature model                   | \~5–10 scalars      |
| Progress, rejection and uncertainty statistics | \~10–20 scalars     |

Start with a deterministic policy for deciding when to extrapolate, correct or explore.

The most revealing test would be an \\(N\\)-dimensional artificial valley with one genuinely soft downhill direction and a strongly curved centerline embedded in many hard coordinates. We can independently control its longitudinal slope, longitudinal curvature, transverse stiffness, turning radius and ambient dimensionality.

If the optimizer can follow that canyon efficiently with one or two evaluations per productive displacement after initial identification, we will have demonstrated the mechanism you are seeking.

We should then deliberately break the assumptions: introduce bifurcations, abrupt changes in direction, multiple soft modes, weakly coupled transverse coordinates, varying stiffness and deceptive local minima.

The really interesting measurement is how many additional evaluations are needed to reorient the optimizer when the valley changes direction.

That separates two qualitatively different abilities: efficient movement along a known valley, and efficient discovery of a new direction when the old one ceases to work.

Finally, I would compare against Powell's NEWUOA/BOBYQA, limited-memory CMA-ES, and simpler gradient-free baselines.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://optimization-online.org\&sz=32)

Optimization Online+1



## My main conclusion

I think the strongest part of your idea is the combination of moving local geometry, predictive extrapolation, and active exploration.

I would not make estimating a general low-rank Hessian the central goal. Instead, try to estimate a much more specific object: a locally curved downhill valley characterized by its tangent, bending direction, longitudinal slope and transverse stiffness.

That directly addresses your original problem of efficiently walking nonlinear canyons.

The optimizer would behave somewhat like a vehicle navigating an unfamiliar road. It remembers the previous direction, estimates how the road bends, accelerates when prediction is reliable, and probes alternative directions when the road suddenly changes. It does not need a global map.

And importantly, this formulation suggests a way to obtain useful information from almost every energy evaluation, including unsuccessful steps.

The next question I would investigate is whether a moving curved-valley model actually predicts future successful steps better than a simpler low-rank quadratic model. That is a concrete and testable distinction. If it does, we have identified an additional source of information—the geometry of the descent trajectory—that conventional momentum heuristics may not be exploiting effectively.


---
---

# USER

OK, so can you write in python or pseudocode some framework we can start with, you can implment several most promising alternativs we wat to try (acording to your judgement) no need to do all the details and souce around it, just the core skelton, and the thing extracted from the paper, so the coding agent which does not have access to these papers have clear starting point to inish the implementation and run the tests/experiments/training...

I can imagine that I will use this branch of FireCore
[https://github.com/ProkopHapala/FireCore/tree/prokop](https://github.com/ProkopHapala/FireCore/tree/prokop)

primarily Electron forcefield (eFF), UFF forcefield, and reactive forcefield RARFF ... as the tcase studies on which I will train it...

---
---

# Chat GPT 6

I created a runnable Python research framework for the gradient-free optimization ideas we discussed, together with a technical handoff for your coding agent. I also inspected the actual `prokop` branch of FireCore to identify the appropriate integration points.

Download ValleyOpt starter package

8 files · Python, NumPy, tests, training and integration notes

## 1. Five experimental optimizers

| Algorithm       | Principle                                                                     | Status             |
| --------------- | ----------------------------------------------------------------------------- | ------------------ |
| `MomentumWalk`  | Remember successful moves; random perpendicular exploration                   | Runnable baseline  |
| `SecantSketch`  | Build an approximate gradient from energy differences, including failed moves | Runnable           |
| `CurvedWalk`    | Estimate and extrapolate the turning direction of the trajectory              | Runnable           |
| `WallProbeWalk` | When hitting a wall, probe transverse directions and estimate slope/curvature | Runnable           |
| `RidgeTR`       | Fit a small quadratic energy model within a moving 2–3D subspace              | Runnable prototype |

All five use energy evaluations only. They have no molecular descriptors, force-field knowledge or neural networks. Their vector-state memory is bounded by \\(O(kN)\\), independently of the trajectory length.

They share a simple `ask()` / `tell()` interface:

```
from valleyopt import RidgeTRopt = RidgeTR(x0, E0, rank=2, step=0.1)for i in range(1000):    x_trial = opt.ask()    E_trial = energy(x_trial)    opt.tell(E_trial)print(opt.best_e, opt.best_x)
```

This is intentionally compatible with a C++ force-field evaluator or future GPU batching.

## 2. Instructions for the coding agent

The most important files are:

Implementation brief

Research objectives, algorithms, extracted ideas from Lion, OMoRF, ASEBO and LM-CMA-ES, numerical caveats, next development tasks and experimental methodology.

Open

FireCore integration

Actual source paths, Python buffer interfaces, eFF example, UFF/MMFF differences, RARFF rotation handling and benchmarking instructions.

Open

The paper notes distinguish established ideas from our own proposals. For example, [OMoRF](https://arxiv.org/abs/2007.04893?utm_source=chatgpt.com) uses separate interpolation sets for discovering a subspace and fitting its reduced quadratic model, while [ASEBO](https://arxiv.org/abs/1903.04268?utm_source=chatgpt.com) adapts exploration within and outside learned active subspaces. These are useful mechanisms to investigate, but the code is not claimed to reproduce either algorithm.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.tandfonline.com\&sz=32)

tandfonline.com+1



## 3. FireCore compatibility

In your branch, I verified that `pyBall/eFF.py` exposes `pDOFs`, `fDOFs` and `eval()`, and that `tests/tEFF/fire_study_lib.py` already contains an instrumented FIRE implementation suitable for comparison.

`pyBall/MMFF.py` also exposes the energy evaluator and coordinate buffers, though MMFF and UFF have different layouts. For RARFF, a dedicated energy-evaluation adapter will be needed.

The package provides `firecore_adapter.py` as the starting point. It handles flattened coordinates, fixed-variable masks, coordinate scaling, and restoring the best configuration.

I have not changed your GitHub repository or run the prototypes against actual FireCore force fields.

## 4. Initial test results

The six unit tests pass. I also ran synthetic benchmarks.

The two-dimensional quadratic-subspace method performs well on Rosenbrock, but the more revealing result concerns dimensionality.

Example: 1,000-dimensional rotated curved valley

180 energy evaluations; initial energy 1.2; global minimum 0

00.000150.00030.00045MomentumSecantCurvedWallRidge

A lower final energy is better; one seed, untuned parameters.

The optimizers barely progress in this high-dimensional example. This is not surprising: discovering an unknown soft direction from a few isotropic scalar probes is extremely difficult.

I also caught and corrected an initially misleading benchmark caused by reusing random seeds for both the landscape and optimizer.

These results are diagnostic, not a meaningful competitive comparison yet.

## 5. What I would develop next

The strongest next experiment is a combination of `RidgeTR` and `WallProbeWalk`:

Keep a 2–3D low-energy subspace and its local quadratic model; when its predictions fail, use carefully selected transverse energy probes to estimate the hard-bath response and discover a better direction.

Then add a tiny controller that decides whether the next evaluation should be spent on extrapolation, a mirror probe, a transverse probe or replacement of a subspace direction.

That should give much more room for algorithmic discovery than simply tuning momentum parameters.

The `train.py` script already provides a basic offline hyperparameter search, so your coding agent can begin with fixed interpretable parameters and only later introduce symbolic program search or a tiny MLP.

To run the package:

```
cd valleyopt
python -m unittest -v
python bench.py --n 12 --neval 1000
python train.py --method ridge --trials 100
```

My suggested first milestone is a reliable benchmark comparing all five candidates against SciPy Powell, a modern derivative-free trust-region solver, and FireCore's native FIRE on eFF. That will tell us which mechanisms are genuinely worth developing before investing in learned control laws.


---
---

# USER

# AGENT



Continuing — let me check whether the methods actually converge with a larger budget, and test a second molecule:\
\
Let me check per-eval wall time (eFF computes forces anyway — an honest caveat from the discussion) and test a perturbed start:\
\
From a strongly perturbed start (E0=-601), \`ridge\` leads at 4.7% gap and reaches 10% in just 79 evals. Let me try the \`load_a\` path — auto-generated electrons, likely a harder initial guess:\
\
Striking result — \`secant\`, \`ridge\`, and \`wall\` all went \*\*below\*\* the FIRE reference (-1616 vs -1259). Either they found a deeper basin or a collapsed-electron pathology. Let me check the final electron sizes:\
\
\`ridge\`'s state isn't converged (|F|=342) but the config looks physical — two tight electrons (\~0.13) on the O core, which eFF wants anyway. Let me verify these are distinct basins by continuing FIRE from ridge's endpoint:\
\
This is a real finding: FIRE from the auto-built start \*\*stalls at -1259\*\* (not converged, just creeping on a plateau), while valleyopt's stochastic descent lands it in the deep basin (-1715). Let me verify on another molecule (CH4):\
\
The auto-generated-electron starts confirm the hypothesis: energy-only methods escape shallow plateaus that stall FIRE. Quick sanity check — is \`secant\`'s -860 on CH4 physical or a collapsed-electron artifact?\
\
The \`secant\` state is physical (tight \~0.17 core electrons, not the 1e-3 floor; FIRE polish continues cleanly to -905.6). I have a complete experimental picture now.\
\
\## Summary\
\
\### Discussion reflection — the fruitful core\
\
The strongest ideas, in order of demonstrated value:\
\
\- \*\*Every evaluation is a gradient projection\*\* (\`SecantSketch\`): even rejected moves yield ΔE = g·s, corrected online in O(N). Cheapest idea, real payoff.\
\- \*\*Curved-valley surrogate\*\* \`E(s,z) = E0 + gs·s + ½ks·s² + ½k⊥(z − Ks²/2)²\` — separating geometric bend K from energy curvatures. Wall penalty ∝ ℓ⁴, so predicting the bend is where the win is. \*\*Not yet implemented\*\* — flagged as the leading next variant.\
\- \*\*Soft subspace + hard bath\*\* (\`RidgeTR\`): 2–3 moving orthonormal directions, O(kN) memory, tiny quadratic fit by ± probes.\
\- \*\*Mirror probes\*\* (\`WallProbeWalk\`): \`E(y±qv)\` separates slope from curvature in 2 evals.\
\- \*\*Direction-search ≠ position-search\*\*: anneal in descent directions while the incumbent stays put.\
\- \*\*Meta-train 5–15 interpretable scalars\*\* by random search on optimization rollouts (not MD trajectories!) — no MLP until fixed coefficients hit their ceiling.\
\
\### Package status — it runs\
\
\- \`python -m unittest\` → 6/6 pass; \`bench.py\` and \`train.py\` work.\
\- \*\*Bug found & fixed\*\*: \`train.py\` sampled \`b/a\` before checking \`a>0\` → \`ZeroDivisionError\` for any method with a zero-lower-bound param (all except ridge). One-line reorder.\
\- Synthetic: \`ridge\` nails 2D problems (Rosenbrock → 9e-25); everything stalls at n≥100 — the k/N information limit, as predicted in the discussion.\
\
\### eFF experiments (\`tests/tEFF/run_valleyopt.py\`, new — runs under ASan build via \`LD_PRELOAD\` + \`CPP_RECOMPILE=0\`)\
\
The adapter works: \`ArrayEnergy(eff.pDOFs, eff.eval, free_mask, positive_indices=esize)\`. Evals cost 14.6 µs — but note eFF computes forces anyway, so energy-only saves nothing on wall time here (a caveat the discussion anticipated).\
\
\*\*On H2O_spins (36 free DOFs, good initial guess):\*\* honest result — Powell dominates (0.64% gap in 2000 evals), FIRE converges, valleyopt methods reach only \~20-27% of the gap in 2000 evals. These prototypes are not competitive in the smooth-descent regime.\
\
\*\*On auto-built electrons (\`load_a\`) — the genuinely interesting result:\*\*\
\
\| mol | start E | FIRE (20k iters) | best valleyopt (3k evals) | + FIRE polish |\
\|---|---|---|---|---|\
\| H2O | -749.7 | -1259 (stalled) | ridge \*\*-1616\*\* | → -1715.3 |\
\| CH4 | -350.9 | -555 (stalled) | secant \*\*-860\*\* | → -905.6 |\
\| CO | -822.4 | -2008 (ok) | secant -2033 | — |\
\
From the bad auto-generated starts, FIRE creeps on a shallow plateau for 20k iterations while the stochastic energy-only methods cross into far deeper basins — verified physical (esize \~0.17-2.8, no floor-hitting), and FIRE-polishable to deep minima. This is the discussion's hypothesis confirmed on real eFF: \*\*the win is in the nonlinear middle of the trajectory, not near the minimum.\*\*\
\
\### Suggested next steps (from IMPLEMENTATION_BRIEF)\
\
1\. Hybrid: valleyopt-style descent → FIRE polish (the -555→-905 CH4 case shows exactly this pipeline working).\
2\. Implement the curved-valley surrogate \`K\`-model — the one major idea not yet coded.\
3\. Soft-direction discovery for RidgeTR when the subspace is exhausted.\
4\. Retrain \`train.py\` parameters on eFF rollouts (the synthetic-tuned params already helped: wall 21.3%→14.9% dE on H2O_spins).\
\
One caveat to keep honest: Powell still beat all prototypes everywhere, and valleyopt stalls badly on synthetic high-D — eFF works because its landscape has exploitable soft collective modes, not because the information limit disappeared.



\---



# USER



few ideas



1\) I do not try to beat BFGS nor FIRE close to minimum. I'm trying to aproach its performance on complex non-convex long-curve-vallye potentials. The stochastic opimized can often jump over small corrugation and barriers, which can make it sometimes faster (especially on non-smooth potentials) than gradient based methods with small number of DOFs.



2\) ad Curved-valley surrogate - yes, I imagine something like that, just I would be carefllu about the qurtic term. Maybe it make sense, but my ideas was that along soft mode we have gradient and quadratic potential, but along the other direction we have still approx quadratic vallyes (resp in our model we estimate curvature in directions perpendicular to our path), these curvature is much stiffer, but still quadratic



3\) Yes the random samples are not suposed to find the gradient in every step. The gradient direction is rather acumulated by the interial history, this is effcient because it amortize previous steps rather than doing many samples each step. But we should also forget old samples (we can forget them by number of iterations which is cheap, but it also make sense to weight them by distance (thrust region). What I was saying about perpendicular directions is something different than tryig to fit exacl curvature, we just do sttistical estimate. Asume there is highest stiffness along unknown direction. If we do  random samples in random direction we can say that stiffness is Kmax\~=K_i/cos(a_i), now imagine we are in low dumensional k-dimensional subspace and there is n-k dimensional stiff sub stapce. We want to only estimate statistically how stiff is this stiff subspace. Basically we do something like low dimenstional (1st,2nd moment) statisticas of many dimenstional hessian. Try to get most of usefull statistical information without doing many samples (and acumulating samples form history), not to find exacl stiffnes along some direction, nor find exact polar basis of hessian. That is important concpeptually, that you I'm for chap and realistc (yet still usefull) goal and not for costly unrealistc goal.&#x20;



4\) we should definitely not vaste to many random samples in wrong direction before we make a move, it is not worth it. We should rather maximally use information form falied samples which we do necessarily anyway. Nice fack is - when we go along good descent direction, most samples are accepted. When the valay and or get kinked in unknond direction, we get stuck and sample randomly around. This is kinda automatic/inherent to this stochastic gradientn-less algorithms, so we should exploit this implicit proiperty of its dynamics. Do not try to do unnecessary failed samples when you descend along correct direction.&#x20;



ad "Direction-search ≠ position-search: when stuck, hold the incumbent fixed and anneal in the space of descent directions — exploration without sacrificing the best configuration."



yes



5\) how to train MPL or even interpretable coefs, what is the loss function for meta-optimization (optimizing the optimizer) ?   Althout total number of iterations is the goal, this give too sparse feedback. Maybe better it to use some proixy e.g. relative drop (|ri|/|ri+k|) after k iteration (sliding window), this is easy to monitor and quntify, or penalize "stuck" situation or unstable oscilation and blow up (mor difficult to qunatify)



\---



reflect on these ideas and update the document that it is clear for future discusstion and development



---

# AGENT



Good clarifications — let me reflect on each, then fold them into \`IMPLEMENTATION_BRIEF.md\` (the living design doc) plus a short results block in \`README.md\`.\
\
\*\*Reflection on your 5 points:\*\*\
\
1\. \*\*Target regime\*\* — agreed, and the eFF runs actually confirmed it: from \`load_a\` auto-generated starts, FIRE stalls on a plateau for 20k iters while stochastic energy-only walks cross barriers into much deeper basins in 3k evals. That's the non-smooth/curved-valley regime, not the parabolic tail.\
\
2\. \*\*Curved-valley surrogate\*\* — right, and it's an important correction to how the doc framed it. All \*fitted\* quantities stay 2nd-order: slope+quadratic along the soft direction, quadratic (much stiffer) \`k⊥\` along transverse directions. The quartic appears only \*implicitly\* through the curved centerline \`z_c(s)=Ks²/2\` — encoded in the moving frame, never fitted as a polynomial term (which would be ill-conditioned).\
\
3\. \*\*Bath statistics, not eigensystem\*\* — this sharpens the goal correctly: for random probe \`v\` in the stiff complement, \`k_v = vᵀHv\` is just a sample from a distribution with mean \`tr(H)/(n−k)\`. If one stiff direction dominates, \`k_i = Kmax·cos²α_i\`, so \`Kmax ≈ max_i k_i/cos²α_i\` is a usable upper-bound statistic even with unknown \`α_i\`. We want 1st/2nd moments accumulated cheaply from history — with forgetting by \*both\* iteration count and distance traveled (samples leave their trust region). Not Hessian reconstruction.\
\
4\. \*\*Exploration rides on organic failures\*\* — the elegant part I'd under-weighted: while descending, most trials are accepted; when the valley kinks, rejections \*automatically\* happen around the incumbent — and those rejected samples ARE the transverse probes, for free. No dedicated probe evals during productive descent.\
\
5\. \*\*Meta-loss\*\* — right, evals-to-convergence is one scalar per rollout (too sparse). Dense proxy: sliding-window relative regret drop \`R_i = regret_i/regret\_{i+w}\` plus secondary stall/oscillation/blowup penalties.\
\
Now the doc updates — surgical edits to \`IMPLEMENTATION_BRIEF.md\`:\
\
Now add the eFF results block to \`README.md\` (including the non-obvious ASan run incantation):\
\
Let me verify the edited sections read correctly:\
\
All edits verified. Here's what I changed:\
\
\## Doc updates\
\
\*\*\`IMPLEMENTATION_BRIEF.md\`\*\* — 4 surgical edits:\
\
1\. \*\*Scientific goal\*\* (top) — added "Target regime (explicit)": we don't aim to beat BFGS/FIRE near minima; target is complex non-convex curved valleys, where stochastic energy-only walks can jump small corrugations/barriers. Referenced the eFF \`load_a\` confirmation.\
\
2\. \*\*Curved-valley surrogate\*\* — your correction is now explicit: \*all fitted quantities stay second-order\* (slope+quadratic along soft dir, quadratic stiffer \`k⊥\` transverse). The \`s⁴\` penalty is implicit in the curved centerline \`z_c(s)=Ks²/2\` — it enters through \*\*frame motion\*\*, never as a fitted polynomial term.\
\
3\. \*\*Priority item 3 → "Hard-bath statistics (moments, not eigensystem)"\*\* — \`k_v = vᵀHv\` from a random probe is just a sample with mean \`tr(H)/(n−k)\`; \`Kmax ≈ max_i k_i/cos²α_i\` as upper-bound statistic. Forgetting by \*both\* iteration count (EMA) and distance traveled (trust-region staleness). Deliberately the cheap-realistic goal vs. Hessian reconstruction.\
\
4\. \*\*Priority item 4 → "Exploration rides on organic failures"\*\* — while descending, most trials are accepted; when the valley kinks, rejections \*automatically\* probe transverse directions for free. Record rejected samples as measurements anchored at the failed point; extra probing only after repeated failures; never burn evals on systematic probes mid-descent.\
\
5\. \*\*Training cost function\*\* — added "Dense feedback — the sparse-signal problem": sliding-window relative drop \`R_i = regret_i/regret\_{i+w}\` as primary signal, plus stuck-fraction and zigzag/blowup penalties, with a starting combined loss form.\
\
\*\*\`README.md\`\*\* — added "eFF observations" section: the results table (H2O/CH4/CO \`load_a\` starts where valleyopt beats stalled FIRE, and the honest Powell-dominates-on-smooth-descent caveat), plus the ASan incantation \`LD_PRELOAD=… CPP_RECOMPILE=0\` needed to run \`tests/tEFF/run_valleyopt.py\` (and a note that Build-opt is stale — missing \`setAutoMass\`).\
\
One open question worth flagging for next discussion: your point 3's \`Kmax ≈ k_i/cos²α_i\` — since \`α_i\` is unknown for a random probe, in practice this means tracking the \*distribution\* of sampled \`k_i\` (max/quantiles) rather than a single estimate. A dominant stiff direction shows up as the upper tail of that distribution; \~O(10) organic samples give a usable estimate of the order of magnitude.


---
---

# Chat GPT 6

This is useful progress. The eFF experiments are particularly interesting because they suggest a more precise target: not faster local relaxation, but faster discovery of productive descent trajectories when gradient-based dynamics becomes trapped in poor basins or narrow nonlinear regions.

I would preserve most of the agent's updated design, but there are three points I would correct before treating it as the basis for further development.

## 1. The hard-bath statistics contain an important mathematical error

The agent wrote:

> \\(K\_{\max}\approx k_i/\cos^2\alpha_i\\) is a usable upper-bound statistic.

This is misleading. Since the angle to the stiffest eigenvector is unknown, we cannot recover \\(K\_{\max}\\) from a measured directional curvature. Indeed, for a positive-semidefinite Hessian,

\\[ k_i=\mathbf v_i^\mathsf T H\mathbf v_i\le K\_{\max}. \\]

Thus the largest observed curvature gives a lower bound, not an upper bound.

But your statistical approach is more interesting than estimating \\(K\_{\max}\\) anyway.

Let \\(H_b\\) be the Hessian restricted to the \\(D=N-k\\) dimensional complementary subspace. For an isotropically random unit vector \\(\mathbf v\\), define

\\[ k_v=\mathbf v^\mathsf T H_b\mathbf v. \\]

Its first two moments are

\\[ \boxed{ \begin{aligned} \mathbb E[k_v]&=\frac{\operatorname{tr}(H_b)}{D},\\\ \operatorname{Var}(k_v)&= \frac{2}{D(D+2)} \left[ \operatorname{tr}(H_b^2) -\frac{\operatorname{tr}(H_b)^2}{D} \right]. \end{aligned} } \\]

This is exactly the kind of cheap information I think you want.

The mean describes the typical stiffness encountered by random transverse exploration. The variance characterizes heterogeneity: is the bath approximately isotropic, or are there exceptional directions with very different stiffness?

You could even estimate an effective spectral rank,

\\[ r\_{\mathrm{eff}}= \frac{(\operatorname{tr}H_b)^2} {\operatorname{tr}(H_b^2)}, \\]

assuming the relevant bath curvatures are nonnegative.

This gives us information about the distribution of stiffness without reconstructing eigenvectors.

One caveat is essential: the formulas assume isotropically sampled directions. Failed proposals are generally not isotropically distributed, because the optimizer selects them based on its current direction and history. Also, a single energy difference includes linear-gradient and soft–hard coupling terms, not just curvature. We must estimate or account for those contributions before interpreting residuals as Hessian measurements.

Therefore, I would make the statistical model of the hard bath, rather than an estimate of maximum stiffness, the central concept.

## 2. The optimization loss should measure anytime performance

I agree with your criticism of evaluations-to-convergence: it is a poor objective for early algorithm discovery, especially on nonconvex problems where different optimizers reach different minima.

However, I would not use the raw sliding-window ratio

\\[ R_i=\frac{E_i-E\_{\rm ref}}{E\_{i+w}-E\_{\rm ref}} \\]

as the primary objective. It becomes unstable near the reference energy and depends strongly on the chosen window. More importantly, optimizing a short-window improvement ratio may incentivize aggressive local descent at the expense of escaping into deeper basins.

I would instead optimize the area under the best-energy trajectory.

Let

\\[ E\_{\rm best}(t)=\min\_{j\le t}E_j \\]

and suppose we have an approximate reference energy \\(E\_{\rm ref}\\) for each training problem.

Define a normalized residual

\\[ r_t= \frac{\max(0,E\_{\rm best}(t)-E\_{\rm ref})} {\max(E_0-E\_{\rm ref},\epsilon)}. \\]

Then use

\\[ \boxed{ L(\theta)= \mathbb E\_{\text{problems}} \left[ \frac1B\sum\_{t=1}^{B} \log(\epsilon+r_t) +\lambda\log(\epsilon+r_B) \right]. } \\]

Here \\(B\\) is the evaluation budget, and \\(\lambda\\) determines the relative importance of final energy versus rapid progress.

This has several advantages:

- Every successful energy improvement contributes to the remaining trajectory's score.
- Early improvements are more valuable than equally large late improvements.
- Temporary failed explorations are not directly punished if they eventually lead to better minima.
- Different potential-energy scales can be normalized.
- No separate artificial oscillation or stagnation penalty is initially necessary.

The reference energy must be estimated carefully, preferably from extensive baseline optimization, and may need updating when a better minimum is found. For synthetic landscapes, we can often know it exactly.

For meta-training fixed coefficients using CMA-ES or random search, the full trajectory still ultimately produces one scalar fitness value. The dense trajectory score improves its informativeness, but doesn't magically provide per-step gradients of the optimizer's parameters.

I would also train over a distribution of evaluation budgets, for example 100, 300, 1000 and 3000 evaluations, to avoid learning an algorithm specialized to only one horizon.

## 3. We should distinguish escaping a bad basin from efficiently following a valley

The eFF results are promising, but I would be careful with the interpretation.

The experiments demonstrate that the stochastic methods find substantially lower-energy configurations than the tested FIRE trajectories from poor initial states.

They do not yet establish that the moving-valley or gradient-estimation mechanism is responsible for this improvement.

A simpler stochastic method might achieve something similar just by making occasional larger jumps.

I would add one particularly important baseline:

FIRE + random perturbations + FIRE, i.e. simple basin hopping.

Use the same total energy/force-evaluation budget. Include random perturbations at several characteristic lengths, and optionally perturb along remembered displacement directions.

If our learned gradient-free optimizer outperforms this simple approach, we have stronger evidence that its accumulation of directional information is genuinely useful.

I would also distinguish two test regimes in the benchmark:

| Regime                                            | What we want to measure                                                   |
| ------------------------------------------------- | ------------------------------------------------------------------------- |
| Smooth curved valley, without competing minima    | Number of evaluations per unit downhill progress; ability to follow bends |
| Corrugated, nonconvex valley with multiple basins | Ability to cross small barriers, escape traps and find deeper minima      |

An optimizer may be excellent at one and poor at the other. We probably want a combination, but understanding each mechanism separately will make the learning problem much easier.

## 4. A further implication for the algorithm architecture

I would now separate the internal state into two components:

\\[ \boxed{ \text{Optimizer state}= \underbrace{\text{trajectory model}}\_{\text{few vectors}} + \underbrace{\text{statistical uncertainty model}}\_{\text{few scalars}} } \\]

The trajectory model remembers successful movement, turning direction and local slope/curvature estimates.

The statistical model tracks how reliably those predictions work, including fluctuations of transverse stiffness and energy-prediction residuals.

The optimizer then has a natural exploration rule:

Exploit while the trajectory model predicts energy changes reliably. Explore when the model ceases to be predictive—not simply when an individual step fails.

This is an important distinction. A single rejected trial could indicate ordinary overshooting, while persistent residuals inconsistent with the local model indicate that the environment has changed.

We can make this decision using explicit statistical rules before introducing an MLP controller.

Overall, I think the development is moving in the right direction. The strongest next experiment would be to determine whether accumulating statistical information from naturally occurring rejected steps improves the prediction of subsequent successful moves.

That is the main mechanism we have hypothesized but not yet isolated experimentally.

I would keep the curved-valley model as a separate experimental variant, rather than making it mandatory. It may turn out that a simpler momentum-based optimizer with a good statistical model of directional uncertainty captures most of the benefit.

And I would correct the \\(K\_{\max}\\) discussion in the living design document: the first and second spectral moments of the hard bath are mathematically well defined, computationally cheap to accumulate under suitable sampling assumptions, and much more closely aligned with your original intention than attempting to infer its stiffest eigenmode.

---
---

# USER

ad "
> \\(K\_{\max}\approx k_i/\cos^2\alpha_i\\) is a usable upper-bound statistic.
This is misleading. Since the angle to the stiffest eigenvector is unknown, we cannot recover \\(K\_{\max}\\) from a measured directional curvature. 
"

I think we should use the assumption of stochasticity! we do not know the angle nor direction, but we can calculate for give n-dimensional sphere what is the statistical distribution of K(angle) and from that estimate best guess. If we have multiple measurements, we can improve that esitmate using bayesian inference.&#x20;

ad "One caveat is essential: the formulas assume isotropically sampled directions. Failed proposals are generally not isotropically distributed, because the optimizer selects them based on its current direction and history. Also, a single energy difference includes linear-gradient and soft–hard coupling terms, not just curvature. We must estimate or account for those contributions before interpreting residuals as Hessian measurements."

example, what I want from you is desing some efficient mathematically sound strategy how to clean these data from the bias. We can take just perpendicular component to the path, and outproject the gradient effect...&#x20;

ad "I would instead optimize the area under the best-energy trajectory."

that is howerver too sparse. I want something which we can train on every small piece of trajectroy (frequent dense signal), not only after trajectroy end, or even ansamble of trajectroies where we exclude some. It is about performance of traning!&#x20;

ad "It becomes unstable near the reference energy and depends strongly on the chosen window. "

The descanti is intitionaly fraction not difference if we assume apporx exponential rate of convergence this fraction should be uniform across the trajectroy.&#x20;

yes, the window should be chosen wisely. Othereise the algorithm can be too greedy and then get stuck way to often. But it will be also penalized for being stuck.&#x20;

But perhaps we should control it by our expected coherence of the descent vectro and memory. We have some memory size resp forgeting rate which caracterize how many steps we think it is worth to keep the information which depend on our assumtion how fast the PES changes. This same asumption tell use what is the relevant number of steps for the window, becase if algorithm is greedy and then stick it means it use the memory inefficiently snd does not employ properly temporal coherence of the underlyain trajectroy resp spatial coherence of the PES.

ad "More importantly optimizing a short-window improvement ratio may incentivize aggressive local descent at the expense of escaping into deeper basins."

Target is not global optimization, if stochastic algorithm jump over barrire, fine, but we should not expect it to do it reliably, local optimizer with grdient like FIRE and BFGS also cannot overcome barrier and escape basins

ad "Exploit while the trajectory model predicts energy changes reliably. Explore when the model ceases to be predictive—not simply when an individual step fails."

reflect on my notes/clarifications, and suggect concrete ideas to implement into the code


---
---

# Chat GPT 6

I agree with your corrections. In particular, I was too dismissive of estimating the maximum stiffness statistically. If we assume a stochastic distribution of probe directions and a reasonably simple prior for the Hessian spectrum, we can infer the likely stiffness of modes that we have never directly identified. We don't need exact reconstruction.

I also agree that the primary training signal should be local in time. We should be able to train the controller from short trajectory segments rather than waiting for complete relaxation.

I see four concrete additions worth implementing:

1. Bayesian hard-bath statistics: infer the distribution of transverse stiffness from randomly oriented trials, not the Hessian eigenvectors.
2. A bias-corrected observation model: separate the expected energy contribution from motion along the known trajectory from the stochastic transverse contribution.
3. Coherence-dependent memory and training horizons: make the relevant history length a property of the locally observed optimization dynamics.
4. Short-horizon counterfactual training: train the controller using comparisons of different decisions starting from the same optimizer state.

The first two can be added to the existing optimizers without changing their basic mechanics. The fourth is what I find most promising for training an MLP or discovering symbolic control rules.

# 1. Bayesian inference of unknown stiff modes

Let's formulate the problem explicitly.

We have an \\(N\\)-dimensional optimization space and a known soft subspace of dimension \\(k\\). The remaining dimension is

\\[ D=N-k. \\]

Suppose the hard-bath Hessian can be approximated by an isotropic background plus a few exceptionally stiff modes:

\\[ \boxed{ H_b=\kappa_0 I+\Delta P_r } \\]

where \\(P_r\\) is the orthogonal projector onto an unknown \\(r\\)-dimensional stiff subspace.

The corresponding eigenvalues are \\(\kappa_0+\Delta\\) in those \\(r\\) directions and \\(\kappa_0\\) elsewhere.

We have no idea where the stiff directions point. But if \\(\mathbf v\\) is uniformly sampled from the \\(D\\)-dimensional unit sphere, the directional curvature is

\\[ q=\mathbf v^\mathsf T H_b\mathbf v =\kappa_0+\Delta z \\]

with

\\[ \boxed{ z\sim\operatorname{Beta}\left(\frac r2,\frac{D-r}{2}\right). } \\]

This follows from the distribution of random projections onto a subspace.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.johndcook.com\&sz=32)

johndcook.com+1



Consequently,

\\[ \begin{aligned} \mathbb E[q]&=\kappa_0+\Delta\frac rD,\\\\[3pt] \operatorname{Var}(q)&= \Delta^2\frac{2r(D-r)}{D^2(D+2)}. \end{aligned} \\]

This is rather powerful: the first two statistical moments of random directional curvature contain information about the strength and number of exceptionally stiff modes, even though we never identify their orientations.

For \\(r=1\\), known \\(\kappa_0\\), and clean samples,

\\[ \hat\Delta\approx D(\overline q-\kappa_0). \\]

The estimate becomes more accurate as samples accumulate. We can do even better by using the entire Beta likelihood instead of just its mean.

For measured curvatures \\(q_1,\ldots,q_m\\),

\\[ p(\kappa_0,\Delta,r\mid q_1,\ldots,q_m) \propto p(\kappa_0,\Delta,r) \prod_i p(q_i\mid\kappa_0,\Delta,r). \\]

This is exactly the Bayesian interpretation you proposed.

We could initially allow only \\(r\in\\{1,2,4,8\\}\\) and use a coarse discretized prior over \\(\kappa_0,\Delta\\). For clean observations, the likelihood is just a transformed Beta density. For noisy ones, convolve it with an error distribution.

The posterior provides not only a best estimate but uncertainty. It can predict, for example, the 90th percentile of transverse stiffness encountered by the next random step.

That predicted percentile is arguably more useful for choosing a safe exploration length than \\(K\_{\max}\\) itself.

### A numerical caveat

With few measurements, \\(\kappa_0\\), \\(\Delta\\) and \\(r\\) may not be separately identifiable. In the limit \\(D\gg r\\), the observations can remain very broadly distributed even when their mean is accurately modeled.

Therefore I would implement two levels:

- Fast runtime estimator: exponentially weighted mean and variance of transverse stiffness, plus observation uncertainty.
- Optional Bayesian estimator: a small posterior over \\(\kappa_0,\Delta,r\\), primarily to test whether spectral assumptions improve predictions.

I would not put numerical integration over a complicated posterior inside every optimization step. The Bayesian model can be updated periodically, or compressed into a few useful predictive quantities.

Also, this spectral model is a hypothesis, not a property of all potential-energy surfaces. We should compare it with a simpler isotropic-bath model and with synthetic Hessians containing broad, continuous spectra.

# 2. Extracting unbiased information from naturally occurring trials

This is the most immediately useful change to the code.

Your suggestion to project perpendicular to the path is right, but there is a subtlety: geometrical projection removes the known soft-subspace displacement, not the energy contribution from an unknown transverse gradient.

We can nevertheless remove that contribution in expectation, using the stochastic symmetry of the exploration directions.

### Derivation

Let \\(U\\) be the current soft-subspace basis. Generate a proposal as

\\[ \mathbf s=U\mathbf a+\delta\mathbf v, \\]

where \\(\mathbf v\\) is sampled uniformly on the unit sphere in the orthogonal complement of \\(U\\).

Locally, expand the energy:

\\[ \begin{aligned} \Delta E \approx{}& \mathbf g_s^\mathsf T\mathbf a+ \frac12\mathbf a^\mathsf TH\_{ss}\mathbf a\\\ &+\delta(\mathbf g_b+H\_{bs}\mathbf a)^\mathsf T\mathbf v\\\ &+\frac12\delta^2\mathbf v^\mathsf TH_b\mathbf v. \end{aligned} \\]

The second line contains the unwanted unknown transverse gradient and coupling between soft and hard modes.

But, conditioned on the optimizer's state and the chosen soft displacement,

\\[ \mathbb E[\mathbf v]=0,\qquad \mathbb E[\mathbf v\mathbf v^\mathsf T]=\frac{I_b}{D}. \\]

Therefore,

\\[ \boxed{ \mathbb E[\Delta E\mid\mathbf a,\delta] = \mathbf g_s^\mathsf T\mathbf a+ \frac12\mathbf a^\mathsf TH\_{ss}\mathbf a+ \frac12\delta^2\mu_b } \\]

where

\\[ \mu_b=\frac{\operatorname{tr}(H_b)}{D}. \\]

This is a very convenient result.

Each ordinary trial becomes a noisy but statistically valid observation of the soft gradient, soft Hessian and average hard-bath stiffness simultaneously.

We do not have to explicitly reconstruct the transverse gradient.

The estimator is accurate to second order in the local expansion; higher-order terms produce finite-step bias.

### Implement this as a tiny online regression

For two soft directions, define the feature vector

\\[ \boldsymbol\phi= \left( a_1,a_2,\frac{a_1^2}{2},a_1a_2, \frac{a_2^2}{2},\frac{\delta^2}{2} \right). \\]

Fit

\\[ \Delta E\approx\boldsymbol\phi^\mathsf T\boldsymbol\theta \\]

with

\\[ \boldsymbol\theta= (g_1,g_2,H\_{11},H\_{12},H\_{22},\mu_b). \\]

This requires just six model coefficients and a \\(6\times6\\) online regression matrix.

I would use recursive least squares, with regularization and a forgetting factor, rather than solve a new least-squares problem from scratch.

The key implementation rule is that every evaluated proposal must enter this estimator, regardless of whether it was accepted or rejected.

### How to avoid selection bias

Three details matter.

First: generate the random transverse component independently of the previous evaluation outcomes, conditional on the optimizer state. Adaptive selection of the amplitude \\(\delta\\) is fine, provided the sampled direction remains symmetric and isotropic for the conditional calculation.

Second: do not update the bath statistics using only failed proposals. Failure selects directions based on their energy response and therefore strongly biases the measured stiffness distribution.

Use accepted and rejected proposals. The fact that the proposal was accepted is an outcome to record, not a filter on the data.

Third: distinguish the stochastic innovation from deterministic steering. If the controller deliberately chooses a particular transverse direction, that displacement is no longer an isotropic sample. Such evaluations are still useful for learning directional structure, but should not enter the isotropic hard-bath estimator as if they were random.

One practical way to enforce this is to construct all proposals from:

\\[ \Delta\mathbf x= \underbrace{U\mathbf a}\_{\text{deterministic exploitation}} + \underbrace{\delta\mathbf v}\_{\text{isotropic stochastic innovation}}. \\]

The optimizer can choose \\(\mathbf a\\) and \\(\delta\\) freely. Only \\(\mathbf v\\) needs the appropriate stochastic symmetry.

### An important distinction about uncertainty

The mean transverse stiffness is identifiable in this setup, but the variance of the one-sided energy residuals has several contributions:

\\[ \operatorname{Var}(\Delta E\mid\mathbf a,\delta) \approx \frac{\delta^2}{D} \\|\mathbf g_b+H\_{bs}\mathbf a\\|^2 + \frac{\delta^4}{4}\operatorname{Var}(q) +\sigma_E^2. \\]

Thus the variance of observed energy residuals is not automatically the variance of directional curvature.

For the first implementation, I would separately track:

- The estimated mean bath stiffness \\(\mu_b\\).
- The prediction-error variance, incorporating unknown transverse gradients and model errors.
- The variance of genuine directional-curvature measurements, when available.

Only the last quantity should feed directly into the simple Beta-based Hessian-spectrum posterior.

### Occasional calibration with mirror probes

Your principle of avoiding unnecessary probes still holds.

But suppose we have been stuck at the same incumbent for several evaluations, and the bath model has large uncertainty.

If we have already evaluated \\(E(\mathbf x+\delta\mathbf v)\\), one additional evaluation at the mirror point gives

\\[ q\_{\mathbf v}\approx \frac{ E(\mathbf x+\delta\mathbf v)+ E(\mathbf x-\delta\mathbf v)-2E(\mathbf x) }{\delta^2}. \\]

The unknown linear gradient cancels. This produces a much cleaner curvature observation, with \\(O(\delta^2)\\) finite-difference error for a smooth potential.

I would trigger these probes only when the model's uncertainty is high and the optimizer has already failed several times. This is consistent with studies of variance reduction using symmetric and orthogonal random finite differences.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://academic.oup.com\&sz=32)

academic.oup.com



# 3. Forgetting by distance, not merely iteration count

I agree that the memory horizon and the training window should originate from the same assumption: the landscape is locally coherent over some characteristic distance.

Let \\(L\\) be the spatial coherence length and \\(\tau\\) a separate iteration-based forgetting timescale.

After an accepted displacement of length \\(\Delta\ell\\), update the accumulated regression statistics using

\\[ \boxed{ \lambda_t= \exp\left( -\frac{1}{\tau} -\frac{\Delta\ell_t}{L} \right). } \\]

For rejected moves, \\(\Delta\ell_t=0\\), since the incumbent has not moved.

This makes old observations decay as the optimizer travels through configuration space.

A few implementation considerations:

- The stored distance should be measured in a consistently scaled coordinate space. Electron widths and atomic Cartesian coordinates, for example, may need different fixed scales.
- Cumulative traveled distance is conservative: it forgets information even if the trajectory returns to an earlier location.
- Regression statistics accumulated in one local frame should be transformed when the frame rotates or the reference point moves. If an accurate transformation is impractical, reset or substantially discount them when the movement exceeds the local trust region.
- The optimizer may use different timescales for momentum and curvature statistics. A useful direction can remain predictive even after the local quadratic approximation becomes inaccurate.

I would make \\(L\\) and \\(\tau\\) trainable symbolic parameters initially.

Later, the controller could adapt them from prediction residuals and successful-step direction correlations.

# 4. Meta-training: your sliding-window approach is worth pursuing

I think your idea of matching the training window to the memory horizon is sensible.

Suppose the optimizer is in a region where its direction and local energy model remain predictive for roughly 20 steps. It makes sense to evaluate its decisions over approximately that timescale, rather than optimize their immediate one-step success or require completion of an entire relaxation.

There is, however, a mathematical issue with using a simple logarithmic ratio.

Let

\\[ r_t=E\_{\rm best}(t)-E\_{\rm ref}+\epsilon S>0 \\]

be the remaining energy gap to a fixed training reference, with \\(S\\) a characteristic energy scale.

If we use

\\[ R_t^{(w)}=\log\frac{r_t}{r\_{t+w}}, \\]

then, when summed over all sliding windows, the logarithmic terms largely telescope. For \\(w=1\\), the total is exactly

\\[ \sum\_{t=0}^{T-1}R_t^{(1)} =\log\frac{r_0}{r_T}. \\]

This means that simply averaging these logarithmic rewards across a whole trajectory does not provide much additional incentive for consistent progress rather than periods of rapid descent followed by stagnation.

Nevertheless, individual short-window rewards remain useful for local policy training.

### A bounded fractional convergence-rate reward

I propose

\\[ \boxed{ R_t^{(w)}= 1- \left(\frac{r\_{t+w}}{r_t}\right)^{1/w}. } \\]

If the remaining gap decreases exponentially,

\\[ r_t=r_0e^{-\gamma t}, \\]

then

\\[ R_t^{(w)}=1-e^{-\gamma}, \\]

independent of both the current progress and the window length.

This has exactly the property you wanted: comparable fractional progress across different phases of the trajectory.

It also has some other useful properties:

- Zero progress gives reward zero.
- A successful decrease gives positive reward.
- An exceptionally large decrease gives diminishing returns rather than an unbounded reward.
- It can be calculated from every short trajectory segment.
- Different horizons can be compared on an approximately per-evaluation basis.

The reference \\(E\_{\rm ref}\\) remains an issue for physical problems. I would use the best-known independently established local reference for each training problem, with a fixed floor \\(\epsilon S\\) to regularize the reward close to convergence. The reference should remain fixed throughout a training segment.

For synthetic problems, we can know the target exactly. For eFF/UFF/RARFF, we can construct a reference database using extensive relaxations. This is an offline training requirement, not information available to the deployed optimizer.

Importantly, this is still local optimization, as you emphasized. We should not reward reliable barrier crossing or global-basin discovery as a separate objective. If crossing a small corrugation improves descent, that is already reflected in the reward.

### Couple the window to the coherence horizon

Let \\(\tau_c\\) estimate the persistence time of useful directional information.

A starting heuristic would be

\\[ w=\operatorname{clip} \left(\lfloor c\\,\tau_c\rceil,w\_{\min},w\_{\max}\right). \\]

For instance, try \\(c=1,2,4\\), with minimum and maximum horizons of perhaps 5 and 100 evaluations.

The coherence estimate might use an exponentially averaged alignment of consecutive successful displacements, combined with local-model prediction errors. We should compare this against simply learning a fixed memory time.

I would also use two horizons during training:

\\[ R_t= 0.7R_t^{(\tau_c)} + 0.3R_t^{(4\tau_c)}. \\]

The exact weights are provisional. The longer horizon discourages excessive greediness while the shorter horizon provides frequent feedback. The horizon lengths should be calculated from the state before evaluating a candidate action, so the controller cannot improve its score by manipulating the evaluation window itself.

## 5. The training procedure I would actually implement: counterfactual trajectory branching

This is where I think we can get a substantial benefit from spending more computation during training without increasing runtime complexity.

Imagine the optimizer reaches a particular state after 100 evaluations. We have its position, incumbent energy, momentum vectors, local quadratic model, bath statistics and coherence estimate.

We can save that state and explore multiple alternative control decisions offline.

For example, try these actions:

| Action A | Long exploitation step along predicted valley |
| -------- | --------------------------------------------- |
| Action B | Short conservative exploitation step          |
| Action C | Moderate stochastic transverse correction     |
| Action D | Large stochastic direction search             |
| Action E | Local-model reset/reorientation               |

Run each alternative for \\(w\\) evaluations, starting from the same saved optimizer state, under the same reference continuation policy after its first decision. Compare their fractional improvements.

This creates an unusually clean supervised-learning signal:

\\[ (\text{dynamical features}) \longrightarrow (\text{which action was best?}). \\]

We can train a tiny MLP, a decision tree, or fit an explicit symbolic rule to these labels.

The training procedure does not require differentiating the force field, differentiating through accepted/rejected moves, or propagating gradients through thousands of optimization iterations.

And although the training process can be expensive, the deployed controller performs only one tiny inference per step.

This is closely related to truncated meta-training and imitation learning approaches studied in learned optimization. Truncation substantially reduces training cost, although short windows can introduce biases that must be tested on complete trajectories.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://learned-optimization.readthedocs.io\&sz=32)

learned_optimization documentation+1



### Why I prefer this to directly optimizing MLP weights by complete rollouts

The controller faces a small number of meaningful decisions, not a high-dimensional action space.

By branching from the same state, we can compare those decisions under nearly identical conditions. This reduces variance in the training labels.

For example, if stochastic direction exploration succeeds from a state where the local model has high prediction error, we can learn exactly which dynamical diagnostics predict that advantage.

We can even inspect whether the controller has learned a simple rule such as:

> Explore only when prediction errors are large, recent attempts have failed, and the transverse uncertainty is substantial.

That is a genuinely interpretable learned heuristic.

We should still validate the resulting policy on complete, unseen trajectories. Counterfactual training optimizes decisions over a limited horizon and can therefore introduce its own myopia.

# 6. Concrete code changes I would give the coding agent

I would implement these in roughly the following order.

## A. Add `OnlineBathModel` to `RidgeTR`

Initially, implement just the six-parameter regression derived above.

The core recursive update is very small:

```
import numpy as npdef bath_features_2d(a, delta2):    """2 soft coordinates + isotropic transverse innovation."""    x, y = a    return np.array([        x, y, 0.5*x*x, x*y, 0.5*y*y, 0.5*delta2    ])class OnlineBathModel:    def __init__(self, prior_var=1e4):        self.theta = np.zeros(6)        self.P = np.eye(6) * prior_var    def update(self, a, delta2, dE, forgetting=0.995):        # Valid only within one consistent local coordinate frame.        phi = bath_features_2d(a, delta2)        pred = phi @ self.theta        self.P /= forgetting        v = self.P @ phi        gain = v / (1.0 + phi @ v)        self.theta += gain * (dE - pred)        self.P -= np.outer(gain, v)        self.P = 0.5 * (self.P + self.P.T)        return dE - pred
```

This is a minimal RLS prototype, not yet a production numerical implementation. Its features must be scaled, its covariance should be regularized, and its state must be reset or transformed when the local frame or reference point changes substantially.

For now, let it provide the soft gradient, soft curvature, mean bath stiffness and predictive uncertainty.

The useful design principle is that the soft model and hard-bath statistics should be estimated jointly, rather than subtracting an imperfect soft model and treating the remainder as exact curvature information.

## B. Add `BathSpectrumPosterior` as a separate experimental module

Implement the model

\\[ H_b=\kappa_0I+\Delta P_r \\]

with a discrete prior over \\(r\\) and simple priors over \\(\kappa_0,\Delta\\).

Initially, feed it clean directional-curvature observations from synthetic tests and occasional mirror probes.

Compare posterior estimates against known eigenvalues.

Only after validating it should we attempt to infer spectral information from one-sided trial-energy residuals using the more complete likelihood that includes unknown transverse-gradient noise.

This separation should prevent us from interpreting random force projections as Hessian heterogeneity.

## C. Add coherence-based forgetting and an exploration controller

The controller can start as a fixed symbolic policy with trainable thresholds.

At each iteration it receives:

```
features = {    "progress_rate":        ...,  # recent fractional improvement    "failure_rate":         ...,  # accepted/rejected history    "direction_coherence":  ...,  # persistence of successful motion    "model_error":          ...,  # normalized prediction residual    "bath_mean":            ...,  # expected transverse stiffness    "bath_uncertainty":     ...,  # prediction uncertainty    "step_length":          ...,}
```

The policy should choose the next proposal mode and step scales.

For example, if the soft model predicts an energy decrease \\(-G\\) and the bath model predicts transverse stiffness \\(q\_{\rm safe}>0\\), we can restrict the exploration amplitude approximately by

\\[ \delta^2\lesssim\frac{2\eta G}{q\_{\rm safe}}. \\]

Here \\(\eta\\) is an allowed ratio of expected transverse penalty to predicted longitudinal gain.

This is a trust-region-style rule, not a hard prohibition on exploration. When the model has become unreliable and the optimizer is stuck, deliberate exploratory moves may need to violate that bound.

Use no extra energy evaluations during ordinary successful descent. Additional probes should become available only in a recovery phase.

## D. Add a short-horizon trainer

I would keep the training code independent of any particular optimizer.

Its conceptual skeleton is:

```
def evaluate_action(snapshot, action, horizon, reference):    opt = restore_optimizer(snapshot)    opt.set_next_action(action)    E0 = opt.best_e    for _ in range(horizon):        x_trial = opt.ask()        opt.tell(opt.energy(x_trial))    scale = max(abs(snapshot.best_e - reference), 1.0)    eps = 1e-6 * scale    r0 = max(E0 - reference, 0.0) + eps    r1 = max(opt.best_e - reference, 0.0) + eps    return 1.0 - (r1 / r0)**(1.0 / horizon)def make_training_example(snapshot, actions, reference):    # Compute horizon from PRE-branch coherence state.    horizon = choose_horizon(snapshot.coherence)    # Clone full optimizer state, not just positions.    scores = [        evaluate_action(snapshot, a, horizon, reference)        for a in actions    ]    return snapshot.scalar_features(), np.asarray(scores)
```

The code is pseudocode for the training interface: `restore_optimizer`, `set_next_action`, `energy`, and `choose_horizon` need implementations. The actual system must also snapshot random-generator state and use matched random streams when comparing actions.

The output need not immediately be a hard label. Training the controller to predict the relative scores of alternative actions may retain more information than simply classifying the winner.

For an initial symbolic optimizer, replace MLP training by directly fitting thresholds or coefficients that maximize the average short-horizon score.

## 7. The experiments that would distinguish the hypotheses

I would run three particularly informative ablations before further architectural expansion.

| Experiment                        | Comparison                                                                                 | What it establishes                                                                                  |
| --------------------------------- | ------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------- |
| Stochastic bath estimation        | True bath statistics vs RLS estimates vs Bayesian estimates                                | Whether we extract useful curvature information from very few evaluations                            |
| Information from natural failures | Momentum only vs momentum + all-sample regression vs regression using rejected trials only | Whether statistically correcting and accumulating the energy observations improves progress          |
| Adaptive exploration              | Fixed exploration amplitude vs uncertainty-controlled amplitude vs tiny learned controller | Whether knowing when to explore is more valuable than optimizing the exploration distribution itself |

For the first experiment, generate synthetic Hessians with known spectra, varying \\(D\\), \\(r\\), background stiffness, gradient magnitude and sample count. Test both clean symmetric curvature probes and ordinary one-sided energy proposals. Plot estimator bias, uncertainty calibration and predictive performance.

For the second, use curved valleys with independently controllable transverse stiffness, bending radius and small corrugations. Keep the evaluation budget identical, and don't let any algorithm access forces.

For the third, use actual eFF, UFF and RARFF configurations, especially states where optimization starts descending productively and subsequently encounters a difficult bend or region of weak progress.

### The main change in philosophy

Previously, we were thinking of the optimizer as something that learns a few good descent vectors.

I think your clarifications suggest a richer but still minimal representation:

\\[ \boxed{ \begin{aligned} \text{Trajectory memory}&:\quad \text{Where does productive motion lead?}\\\ \text{Local energy model}&:\quad \text{What improvement should we expect?}\\\ \text{Bath statistics}&:\quad \text{How risky is unexplored transverse motion?}\\\ \text{Coherence memory}&:\quad \text{How long should we trust what we learned?} \end{aligned} } \\]

These are four different kinds of information, but none requires learning molecular structure or storing a large history.

My priority would be the joint soft-subspace/hard-bath regression first, followed by counterfactual short-window training. The former extracts more information from evaluations we already pay for; the latter gives us an efficient way to discover whether and when that information improves decisions.

The Bayesian spectral model is mathematically appealing and worth implementing, but I would regard it as an experimental extension until we demonstrate that predicting the distribution of transverse stiffness actually leads to better optimization trajectories—not merely better estimates of Hessian statistics.

---
---

# USER

OK, we tried the gradient-less method and it is not bad but Powell and FIRE is much faster so not much gain...

One option would be try to beat Powell by imporving it.

But perhaps more frutful is to focus on imporving FIRE, by a bit more complex learned heuristics...

---

from my coding agent

Let me check the topical audit structure first, then clarify metrics and make plots:

Now the metrics question — you're right, a single `dE%` number doesn't show convergence *speed*. `dE% = 100·(E_best − E_ref)/(E0 − E_ref)` — the fraction of the start→reference gap still remaining at budget end. To show speed we need the full best-energy-vs-evals curve. Let me add trace saving + plotting to the benchmark:

NPZ saved but no PNG — plotting must have failed. Check:

ASan's `LD_PRELOAD` breaks matplotlib's `ft2font` C extension (missing `__cxa_throw` interceptor for non-instrumented libs). Fix: spawn plotting in a child process with `LD_PRELOAD` stripped:

The module-level `pyBall` import loads the ASan lib even for `--plot`. Make the imports lazy:

Plots generated. Let me view one to check it conveys convergence speed:

The data confirms the plots show convergence speed correctly. Check where to keep the PNGs (they're gitignored):

PNGs are gitignored runtime artifacts — keep them in tests/tMMFF and document regeneration. Now write the comprehensive audit report:

Now cross-reference from the relevant audit docs:

Add the convergence plots reference to the valleyopt README:

Done — everything is written up and plots generated.

## The "% of what" — clarified

`dE%` = `100·(E_best − Eref)/(E0 − Eref)` — the fraction of the start→reference gap *still remaining* at budget end. It's a single end-point number — you're right that it says nothing about *speed*. Two things show speed:

- **Milestone columns** (`0.5 0.1 0.02 0.005`) = first eval index reaching that relative level — evals-to-threshold.
- **Convergence curves** — now generated: `run_valleyopt.py` saves `<mol>_convergence.npz` + `.png` in `tests/tMMFF/` (regenerate with `--plot <name>`; plotting runs in a clean child because ASan's `LD_PRELOAD` crashes matplotlib's ft2font).

## What the curves show (xylitol, gap remaining vs evals)

| method | @50 | @200 | @500 | @2000 |
|---|---:|---:|---:|---:|
| FIRE* (forces) | 0.69 | 0.33 | 0.21 | **0.044** |
| powell | 0.98 | 0.79 | 0.58 | 0.240 |
| secant | 0.91 | 0.73 | 0.52 | **0.306** |
| ridge | 0.92 | 0.77 | 0.63 | 0.376 |
| curved/wall/momentum | ~0.99 | ~0.9 | ~0.8 | 0.54–0.72 |

Reading: FIRE is ~4-5× more eval-efficient everywhere. `secant`'s amortized sketch descends steadily (log-linear-ish); `powell` starts slower but has the steepest tail; the random-direction prototypes (momentum/wall/curved) essentially stall at 66 DOFs. Same story at 39/168 DOFs.

## Audit report

New file: **`doc/topical_audit/MMFF_energy_eval_debug_2026-11-08.md`** — comprehensive: DOF layout reference (`ff` vs `ffl` layouts, dead pi DOFs), all 7 bugs with symptoms/fixes/locations, the GPU precedent (`updateAtomsMMFFf4` pi normalization), remaining caveats (PBC clamp divergence, repeated-init UAF, ASan+matplotlib, UFF ndims unchecked), metric definitions, and the verification-gate command. Cross-linked from `topical_audit.md` §2a, `intramolecular_forcefields.md` §2 ("READ FIRST" block), and `nonbonding_forcefields.md` (half-share partition convention).

---
---

# Chat GPT 6

I agree with changing direction. The gradient-free experiments have been useful because they showed where the information bottleneck lies. With FIRE, we already obtain an entire force vector at every evaluation, and the challenge is to use the force history more intelligently without paying for an expensive Hessian approximation.

I think there is room for improvement beyond simply learning FIRE's timestep and damping coefficients. In particular, I would investigate three mechanisms:

1. Adaptive inertia and preconditioning: allow different dynamical modes to relax at different rates.
2. Predictive steering: use the change in force to anticipate bends in a valley instead of waiting for FIRE's velocity–force alignment to deteriorate.
3. A learned dynamical controller: recognize persistent downhill motion, transverse oscillation, overshooting and changing curvature, then choose how aggressively to continue.

One useful finding from inspecting your `prokop` branch of FireCore: you already have some relevant experiments in `tests/tEFF/fire_study_lib.py`, including `fire_groups()`, `fire_secant()`, mass preconditioning and a diagonal-Newton reference. `DynamicOpt.cpp` also contains `autoMass`, `massMetric`, and limited-memory gradient-history methods.

We should build on these rather than repeat them. In fact, the prior mass-adaptation experiments are potentially more important than adding an MLP immediately, because they test whether modifying the dynamics in different directions is beneficial at all.

## 1. First: there are already improved FIRE algorithms worth testing

I found three particularly relevant developments.

| Algorithm                                                                                 | Main modification                                                                               | Relevance                                                      |
| ----------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| [FIRE 2.0](https://doi.org/10.1016/j.commatsci.2020.109584?utm_source=chatgpt.com) (2020) | Improved integration, half-step-back on overshoot, safer timestep control                       | Essential reference baseline                                   |
| [ABC-FIRE](https://doi.org/10.1016/j.commatsci.2022.111978?utm_source=chatgpt.com) (2023) | Bias-corrected velocity mixing, reducing slow restarts                                          | Cheap modification with reported iteration reductions of 9–65% |
| [AARE](https://arxiv.org/abs/2410.21837?utm_source=chatgpt.com) (2024)                    | Avoids unnecessarily zeroing velocity; uses CG-like directions and angle-dependent acceleration | Particularly close to your desired approach                    |

The AARE paper is especially interesting. Instead of completely killing velocity whenever FIRE detects overshooting, it changes the update rule depending on the angle between the force and the previous direction. It switches among Polak–Ribière, Fletcher–Reeves, Hestenes–Stiefel and steepest-descent directions. The reported benchmarks include improvements exceeding fivefold in some small-dimensional examples, though this doesn't establish similar improvements for eFF or MMFF.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arxiv.org+1



Another relevant paper is [Enhanced FIRE for multiscale simulations](https://doi.org/10.1016/j.commatsci.2024.113234?utm_source=chatgpt.com) (2024). It explicitly addresses different relaxation rates across stiff and soft degrees of freedom by choosing effective masses.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com\&sz=32)

sciencedirect.com



These papers suggest that improving FIRE is not necessarily a matter of finding better global coefficients. Changing how momentum responds to different dynamical modes may be more fruitful.

## 2. My preferred experiment: selective quenching and low-rank adaptive inertia

Ordinary FIRE treats the velocity as a single vector. When its power

\\[ P=\mathbf v\cdot\mathbf F \\]

becomes negative, the velocity is quenched.

But imagine a narrow valley where most of the motion is productive along the valley, while one stiff transverse mode oscillates.

The transverse oscillation can eventually make \\(P\\) negative, causing FIRE to discard both the unwanted transverse velocity and the useful longitudinal velocity.

That is an obvious opportunity for a slightly smarter optimizer.

Let

\\[ \mathbf u=\frac{\mathbf s}{\\|\mathbf s\\|} \\]

be the recent successful direction.

Decompose velocity and force into components parallel and perpendicular to \\(\mathbf u\\):

\\[ \begin{aligned} \mathbf v&=\mathbf v\_\parallel+\mathbf v\_\perp,\\\ \mathbf F&=\mathbf F\_\parallel+\mathbf F\_\perp. \end{aligned} \\]

Now we can adapt these components independently.

For example, if

\\[ \mathbf v\_\perp\cdot\mathbf F\_\perp<0 \\]

while

\\[ \mathbf v\_\parallel\cdot\mathbf F\_\parallel>0, \\]

we could strongly damp the transverse velocity without discarding longitudinal momentum.

This is much closer to your idea of preserving useful inertia while suppressing oscillations in stiff modes.

But we can do even better than splitting into just parallel and perpendicular components.

### Use force changes to identify an important transverse mode

From two successive force evaluations,

\\[ \Delta\mathbf F=\mathbf F_t-\mathbf F\_{t-1}, \qquad \mathbf s=\mathbf x_t-\mathbf x\_{t-1}. \\]

We know approximately

\\[ \Delta\mathbf F\simeq-H\mathbf s. \\]

The force change perpendicular to the displacement,

\\[ \mathbf q= \Delta\mathbf F- \frac{\Delta\mathbf F\cdot\mathbf s} {\\|\mathbf s\\|^2}\mathbf s, \\]

provides a direction in which the Hessian couples to our motion.

This is not necessarily an eigenvector of the stiffest mode. Nevertheless, it can identify a direction responsible for a significant force response.

We can normalize \\(\mathbf q\\) and use it to define one dynamically discovered mode.

Then apply selective damping or an effective mass along that mode.

For example, a simple positive-definite inverse-mass operator would be

\\[ \boxed{ M^{-1}=I-(1-\mu)\hat{\mathbf q}\hat{\mathbf q}^{\mathsf T}, \qquad 0<\mu\le1. } \\]

When \\(\mu\\) is small, acceleration along \\(\mathbf q\\) is suppressed, while acceleration in orthogonal directions remains unchanged.

We could keep two or three such orthogonalized directions, replacing them as the trajectory evolves.

This is a low-rank dynamical preconditioner, requiring only a few full-size vectors, dot products and vector additions.

Unlike a diagonal mass model, it is rotation-equivariant and does not require knowing which coordinates belong to atoms, electron widths or other physical variables.

The actual stiffness and relevance of each direction must be estimated from history, and the mode should be forgotten when it becomes stale.

## 3. Second experiment: predictive FIRE using force differences

Another possibility is to predict the force at the next position.

For approximately smooth motion,

\\[ \mathbf F\_{t+1}\approx \mathbf F_t+ \beta_t(\mathbf F_t-\mathbf F\_{t-1}), \\]

where \\(\beta_t\\) accounts for the anticipated displacement relative to the previous one.

Instead of using only the present force for steering, use a limited extrapolation:

\\[ \mathbf F\_{\rm steer}= \mathbf F_t+\lambda_t\Delta\mathbf F. \\]

This resembles momentum-based prediction, but applied to the force field rather than just the optimizer velocity.

The controller can learn whether extrapolation is useful based on the previous prediction errors.

When the force evolves smoothly, extrapolation may help anticipate a valley turn. When forces oscillate rapidly, it should be suppressed.

There is an obvious danger: extrapolating a stiff oscillatory force may amplify rather than damp instability. Therefore I would constrain \\(\lambda_t\\), measure prediction accuracy, and keep the actual force for descent and safety checks.

This can be tested with no additional force evaluations.

## 4. Third experiment: learned controller for FIRE's internal parameters

Only after testing the previous mechanisms would I introduce an MLP.

Its inputs should be a few dimensionless dynamical observables:

| Input                                                        | Interpretation                             |
| ------------------------------------------------------------ | ------------------------------------------ |
| \\(\cos(\mathbf v,\mathbf F)\\)                              | Current downhill alignment                 |
| \\(\cos(\mathbf F_t,\mathbf F\_{t-1})\\)                     | Force-direction persistence                |
| \\(\cos(\mathbf s_t,\mathbf s\_{t-1})\\)                     | Trajectory coherence                       |
| \\(-\mathbf s\cdot\Delta\mathbf F/\\\|\mathbf s\\\|^2\\)     | Directional stiffness, after normalization |
| \\(\\\|\Delta\mathbf F\_\perp\\\|/\\\|\Delta\mathbf F\\\|\\) | Strength of transverse force coupling      |
| Recent velocity-quench frequency                             | Persistent overshooting                    |
| Predicted versus actual force change                         | Reliability of the local model             |

The controller would output a small number of coefficients, such as:

\\[ (\alpha,\Delta t,\mu,\lambda,\gamma\_{\rm quench}). \\]

Here \\(\mu\\) controls the effective mass of the identified mode, \\(\lambda\\) controls force extrapolation, and \\(\gamma\_{\rm quench}\\) determines how much velocity to retain during an overshoot.

I would not allow the controller to propose arbitrary directions. It should select coefficients of the explicit update rule.

And I would initially train only the coefficients of a symbolic policy before introducing the MLP.

## 5. One immediate concern with the FireCore implementation

From the current `DynamicOpt.cpp`, the main `move_FIRE()` implementation still uses conventional velocity mixing followed by a leapfrog-style step. I did not see FIRE 2.0's explicit half-step-back correction on negative power.

There is also an existing mass estimator using

\\[ k_b\approx\frac{\\|\Delta\mathbf F_b\\|}{\\|\Delta\mathbf x_b\\|}, \\]

followed by

\\[ M_b^{-1}=\frac{k\_{\rm ref}}{k\_{\rm ref}+k_b}. \\]

That is a sensible cheap heuristic, but it measures a response involving coupling to other coordinates, not necessarily the intrinsic stiffness of block \\(b\\).

Therefore, the most useful immediate experiments would be:

| Variant                         | Scientific question                                                                |
| ------------------------------- | ---------------------------------------------------------------------------------- |
| FIRE 2.0 / ABC-FIRE             | How much improvement is available simply by refining existing FIRE dynamics?       |
| FIRE + existing `autoMass`      | Does blockwise adaptive inertia help on the selected test problems?                |
| FIRE + selective quenching      | Can we preserve productive longitudinal momentum during transverse oscillations?   |
| FIRE + rank-1/2 mass adaptation | Can we suppress stiff collective modes without problem-specific coordinate groups? |
| FIRE + predictive steering      | Does extrapolating force changes improve navigation of curved valleys?             |
| Learned controller              | Can we combine these mechanisms more effectively than fixed thresholds?            |

I would measure not only energy versus force evaluations but also the number of velocity resets, time spent at the minimum timestep, directional coherence, overshoot frequency, and rate of decay of the force norm.

The last point is important: a force-based optimizer can reduce energy quickly while remaining far from a stationary configuration.

For eFF, where the agent measured approximately 14.6 µs per evaluation, the prototype should eventually run inside C++, not as a Python loop with multiple array allocations and small MLP evaluations. Otherwise, optimization overhead may dominate the very force-field calculation we are trying to accelerate.

## 6. Where I would place the main research bet

First priority: mode-selective FIRE

My strongest candidate. Separate productive momentum from stiff oscillatory motion using one or two dynamically estimated modes. Test simple, interpretable damping and effective-mass rules before learning them.

Second priority: predictive force steering

Test whether extrapolating force changes helps on strongly curved synthetic valleys. Measure its reliability and combine with selective damping only if independently useful.

Third priority: learned control law

Train a controller on short trajectory windows, using the same fractional-progress objective we previously discussed. Optimize controller coefficients or symbolic expressions first, and compare against a tiny MLP later.

For training, I would retain your notion of a coherence-dependent window. With force information, we can additionally measure whether the force norm decreases over that window, without introducing extra evaluations.

The deeper idea is that FIRE already knows how to accelerate down a valley, but its response to stiff transverse motion and changing valley geometry is relatively crude. Adding a little memory of the force response rather than only the velocity could be valuable.

The main uncertainty is whether the identified modes remain coherent long enough to justify their extra state and control parameters. A simple test on synthetic curved valleys with known stiffness ratios, followed by eFF/UFF/RARFF comparisons, should answer that before we invest in learned controllers.

I would be more optimistic about this direction than gradient-free optimization, especially because all the additional information can be extracted from forces FireCore already computes.

---
---

# USER

OK, I have to say we should not imporve only on FIRE but also on Anderson, or at least make cheaper simple version comparable to anderson without storing large history of vectors...

One odea could be to have per-atom mass shich is adjusted based on measurements along the trajectroy (that can be permanent nor forgeting) which can replace variable timestep used curencly in FIRE

read this document reflect on it, and suggest what we should try

---
---

# Chat GPT 6

After reading the full `eFF_optimizer_speedup.md`, especially §§6–8, I think we have a clearer direction than simply modifying FIRE.

I would try to develop a family of optimizers between FIRE and Anderson: retain FIRE's cheap inertial dynamics, add persistent adaptive stiffness information, and incorporate one or two secant corrections that capture some of Anderson's acceleration without storing a long history.

Your idea of dynamically adjusted per-atom masses fits this approach. But the experiments in the document also tell us where mass adaptation alone is insufficient.

## 1. What the experiments actually establish

The most relevant results are in §8.9, after the Anderson history bug was corrected. These supersede several earlier tables in the document.

Corrected frozen-core H₂O benchmark

Mean force evaluations across original, rotated and distorted configurations; raw free-DOF force tolerance \\(10^{-3}\\)

| Method                                        | Free nuclei | Fixed nuclei |
| --------------------------------------------- | ----------- | ------------ |
| Tuned FIRE                                    | 184         | 67           |
| FIRE + per-particle secant masses             | 137         | 89           |
| BB1 spectral gradient                         | 662         | 45           |
| Anderson m=5–8                                | 349         | 33           |
| Preconditioned FIRE (60 steps) → Anderson m=8 | 108         | —            |
| Scaled dense BFGS                             | \~72        | —            |

eFF_optimizer_speedup.md

Three conclusions are particularly important.

First, adaptive particle masses already work. The tested scheme uses one stiffness estimate per particle's three Cartesian coordinates, and a separate stiffness for each electron width. It provides a meaningful improvement in the free-nuclei case and is naturally rotation-invariant. However, it actually makes the fixed-nuclei case slower.

Second, Anderson's advantage is not simply a better scalar step size. BB1 already estimates a characteristic curvature from consecutive gradients, but performs extremely poorly on the coupled free-nuclei system. Anderson exploits correlations between several displacement and gradient-change vectors. This is information that diagonal masses cannot fully reproduce.

Third, the combination performs better than either method separately. Anderson works poorly from the initial free-nuclei configuration, but becomes effective after preconditioned FIRE approaches the relevant basin.

For me, the third result is the strongest motivation for a new algorithm. Instead of explicitly switching between two separate optimizers, perhaps we can create a simple dynamical rule that gradually acquires Anderson-like behavior as its local secant information becomes reliable.

One qualification: these experiments are still mostly variants of one small H₂O system. We have evidence for a useful mechanism, not yet for general superiority.

## 2. Can adaptive per-atom masses replace FIRE's global timestep?

Largely yes, but I would formulate it as adaptive mobility rather than physical mass. This removes an unnecessary distinction between timestep and mass for a relaxation algorithm.

Consider simple inertial dynamics:

\\[ \mathbf v\_{t+1}=\beta\mathbf v_t+\Delta t\\,M^{-1}\mathbf F_t, \qquad \mathbf x\_{t+1}=\mathbf x_t+\Delta t\\,\mathbf v\_{t+1}. \\]

Introduce the displacement memory \\(\mathbf p=\Delta t\\,\mathbf v\\). For constant \\(\Delta t\\),

\\[ \boxed{ \begin{aligned} \mathbf p\_{t+1}&=\beta\mathbf p_t+B_t\mathbf F_t,\\\ \mathbf x\_{t+1}&=\mathbf x_t+\mathbf p\_{t+1}, \end{aligned} \qquad B_t=\Delta t^2 M_t^{-1}. } \\]

This makes the important point explicit: for numerical optimization, the effective mobility \\(B\\), not mass and timestep separately, determines the force-driven displacement.

We could dispense with an explicit time variable and adapt one scalar \\(b_i\\) per degree of freedom.

FIRE's normalized velocity–force mixing and power-based quenching can still be applied to \\(\mathbf p\\).

### What should we learn for each particle?

I would start with the already-tested trajectory secant:

\\[ \kappa_b^{\rm obs}= \frac{\\|\Delta\mathbf F_b\\|} {\\|\Delta\mathbf x_b\\|+\epsilon}. \\]

Here \\(b\\) denotes a three-dimensional particle block, not an individual Cartesian component.

Then define the mobility

\\[ \boxed{ b_b=\frac{b_0}{1+(\kappa_b/\kappa\_{\rm ref})^p}. } \\]

The exponent \\(p=1\\) reproduces the general behavior of the existing adaptive-mass implementation. Testing \\(p=0.5,1,1.5\\) would tell us whether full inverse-stiffness scaling is too aggressive.

For eFF, keep electron widths as separate scalar blocks.

The more interesting choice is how to accumulate the stiffness estimate.

### Persistent memory versus forgetting

I would test three options:

| Memory rule                   | Advantage                                             | Potential problem                                               |
| ----------------------------- | ----------------------------------------------------- | --------------------------------------------------------------- |
| Ordinary EMA                  | Adapts continuously to changing curvature             | Can forget an important stiff mode too quickly                  |
| Maximum encountered stiffness | Retains a conservative stability estimate             | One anomalous observation can permanently immobilize a particle |
| Fast increase, slow decrease  | Remembers stiff modes but gradually recovers mobility | Introduces two timescales                                       |

The third is my preferred candidate.

For example, work with logarithmic stiffness:

\\[ \begin{aligned} z_b&=\log(\kappa_b+\epsilon),\\\ z_b&\leftarrow z_b+ \eta\_\pm\left(\log(\kappa_b^{\rm obs}+\epsilon)-z_b\right), \end{aligned} \\]

with \\(\eta\_+=0.2\\) when the observation indicates greater stiffness and \\(\eta\_-=0.01\\) otherwise.

Those numbers are experimental starting values.

This creates a useful asymmetry: detecting unexpectedly stiff dynamics rapidly increases mass, but an isolated soft observation does not immediately remove that protection.

I would also test learning masses during an initial calibration phase, then freezing them. For approximately stationary problems, fixed learned masses may be preferable to permanently modifying the dynamics.

The trajectory secant is still not a true diagonal Hessian estimate: coupled motion can change a particle's force even when that particle barely moves. To avoid spurious huge stiffness estimates, reject observations with insignificant local displacement and test an EMA ratio of accumulated \\(\\|\Delta F_b\\|^2\\) and \\(\\|\Delta x_b\\|^2\\), rather than always dividing each individual pair.

### Why this cannot completely eliminate global step control

For a harmonic potential with Hessian \\(H\\), the stability of the undamped inertial update depends on the eigenvalues of

\\[ B^{1/2}HB^{1/2}. \\]

For the simple leapfrog discretization, stability requires its largest eigenvalue to be below 4.

A diagonal mobility may reduce this spectral radius dramatically, but it cannot generally eliminate coupled stiff modes. Even if all diagonal elements of \\(H\\) are moderate, a collective mode may still be very stiff.

Therefore I would use constant nominal mobility during productive descent, but retain an emergency global scaling or trust-radius safeguard when instability is detected.

The distinction from ordinary FIRE is that the global multiplier would no longer have to increase and decrease continually to compensate for heterogeneous stiffness. Most adaptation would occur in the individual mobilities.

## 3. The most promising cheap alternative to Anderson: one-secant correction

This is the experiment I would prioritize.

Anderson currently stores 5–8 pairs of displacement and gradient-change vectors. But Anderson with history length one has a particularly simple explicit form, closely related to a rank-one inverse-Jacobian secant update. The relationship between Anderson acceleration and multisecant/Broyden methods is well established.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://epubs.siam.org\&sz=32)

epubs.siam.org+1



Define

\\[ \mathbf g_t=-\mathbf F_t,\qquad \mathbf s=\mathbf x_t-\mathbf x\_{t-1},\qquad \mathbf y=\mathbf g_t-\mathbf g\_{t-1}. \\]

Let the baseline approximate inverse Hessian be \\(\eta D\\), where \\(D\\) is our diagonal adaptive mobility/preconditioner.

Then propose

\\[ \boxed{ \begin{aligned} \gamma&=\frac{\mathbf y\cdot\mathbf g_t} {\mathbf y\cdot\mathbf y+\lambda},\\\\[2pt] \mathbf d&= -\eta D\mathbf g_t -(\mathbf s-\eta D\mathbf y)\gamma. \end{aligned} } \\]

This is a damped one-secant Anderson-like update.

Every term has a clear interpretation:

- \\(-\eta D\mathbf g\\) is ordinary preconditioned gradient descent.
- \\(\mathbf s-\eta D\mathbf y\\) measures the discrepancy between actual motion and what the baseline mobility would predict from the force change.
- \\(\gamma\\) determines how strongly that historical correction is relevant to the current gradient.
- \\(\lambda\\) regularizes poorly determined corrections.

With \\(\lambda=0\\), the corresponding rank-one inverse-Jacobian approximation satisfies the secant equation exactly, provided \\(\mathbf y\ne0\\).

The evaluation cost is essentially zero beyond the force evaluation already performed. The update needs a few dot products and vector additions, with no matrix inversion.

Importantly, this is not a new mathematical discovery. It is a minimal rank-one multisecant rule. The research question is whether combining it with persistent adaptive mobility and FIRE-like inertia produces a useful optimizer.

### Why this could work better than ordinary BB1

BB1 estimates one characteristic inverse curvature:

\\[ \eta\_{\rm BB1}= \frac{\mathbf s\cdot\mathbf s} {\mathbf s\cdot\mathbf y}. \\]

It can change the magnitude of a gradient step, but cannot change its direction.

The rank-one correction above can modify the direction of the proposed displacement according to historical force response. This gives it a way to compensate for coupled modes.

But one secant direction cannot represent several independently stiff modes. I would therefore compare:

- Rank 1: one recent secant pair.
- Rank 2: two secant pairs with a \\(2\times2\\) regularized solve.
- Rank 3: three pairs.
- Existing Anderson m=5 and m=8.

The rank-2 version is still extremely cheap. We only need four history vectors and a tiny system solve.

Secant-history memory for one million coordinates

Float64, displacement and gradient-difference vectors only

0 MB35 MB70 MB105 MB140 MBRank 1Rank 2Rank 3Anderson m=5Anderson m=8

Excludes current positions, forces and other optimizer state.

The important experiment is whether m=2 captures most of m=5–8's benefit. If it does, there is little reason to retain the larger history.

## 4. A potentially better combination: inertia plus a secant correction

Rather than switching between FIRE and Anderson after 60 evaluations, we could attempt to continuously combine their update mechanisms.

Let \\(\mathbf p\_{\rm FIRE}\\) be the displacement proposed by an adaptive-mobility FIRE step, and \\(\mathbf d\_{\rm secant}\\) the rank-one or rank-two correction described above.

Then

\\[ \boxed{ \Delta\mathbf x= (1-w_t)\mathbf p\_{\rm FIRE} +w_t\widetilde{\mathbf d}\_{\rm secant}, } \\]

where \\(\widetilde{\mathbf d}\_{\rm secant}\\) is scaled or capped to a compatible displacement length.

The central question becomes how to determine \\(w_t\\).

I would start with an explicit reliability measure, not an MLP.

### Can we predict how useful a secant correction will be?

One simple idea is to compare the observed gradient difference against the response predicted by the existing local approximation.

If consecutive secant corrections become coherent and the predicted energy reduction agrees with the actual reduction, increase \\(w_t\\). If the model becomes unreliable or proposes uphill directions, decrease it.

For example, calculate

\\[ c\_{\rm descent}= -\frac{\mathbf g_t\cdot\mathbf d\_{\rm secant}} {\\|\mathbf g_t\\|\\,\\|\mathbf d\_{\rm secant}\\|+\epsilon}. \\]

If this becomes negative, the secant direction is locally uphill, so we should reject or substantially suppress it.

But positive alignment alone is insufficient. A direction can be downhill while having a disastrously large step.

We should therefore also record the predicted versus actual energy improvement and use a trust-radius restriction.

My preferred implementation would initially give the secant correction a small weight, increasing it only after several reliable steps. On a discontinuity or severe prediction failure, reset its history without necessarily resetting the FIRE momentum.

This could avoid the poor early behavior of Anderson while gradually acquiring its acceleration near smoother regions.

There is a risk that this interpolation combines the weaknesses of both optimizers rather than their advantages. That is why I would test pure rank-1/2 Anderson first, then the blended variant.

## 5. Another cheap possibility: identify one stiff collective mode

The study reveals an important limitation of per-particle masses.

In all-electron eFF, the stiffest modes can involve coordinated motion of a nucleus and its bound core electron. Giving both particles large masses slows the relative stiff motion, but also slows their common translation.

This is an inherent limitation of diagonal preconditioning.

A rank-one preconditioner can distinguish these two kinds of motion without knowing anything about molecular geometry.

We already measure the force difference

\\[ \mathbf y=-\Delta\mathbf F\approx H\mathbf s. \\]

If a strongly excited stiff mode dominates the force change, the direction of \\(\mathbf y\\) often approximates that mode.

We could track its normalized direction \\(\mathbf u\\), using a smoothed history of reliable force changes, and define

\\[ \boxed{ B=D^{1/2} \left(I-\rho\mathbf u\mathbf u^\mathsf T\right) D^{1/2}, \qquad 0\le\rho<1. } \\]

Here \\(D\\) is the diagonal adaptive mobility.

This suppresses acceleration along one empirically discovered stiff direction while preserving the other directions.

Applying \\(B\\) to a force requires just one extra dot product and a vector correction.

The same idea could be extended to two orthogonal modes.

It is particularly attractive for eFF's stiff coupled electron–nucleus motion, although I would not assume the force-difference direction always identifies the largest-curvature eigenvector. It may instead identify the mode most strongly excited by the current trajectory.

A good test would compare this rank-one preconditioner with ordinary per-particle mass adaptation on a deliberately rotated, coupled quadratic potential. If it cannot outperform diagonal masses on such a controlled example, I would not pursue it further.

## 6. A scalar diagnostic worth extracting: BB1 versus BB2

We have another potentially useful source of information that requires no extra state.

Calculate

\\[ \eta_1=\frac{\mathbf s^2}{\mathbf s\cdot\mathbf y}, \qquad \eta_2=\frac{\mathbf s\cdot\mathbf y}{\mathbf y^2}. \\]

For positive \\(\mathbf s\cdot\mathbf y\\),

\\[ \boxed{ \frac{\eta_1}{\eta_2} = \frac{\\|\mathbf s\\|^2\\|\mathbf y\\|^2} {(\mathbf s\cdot\mathbf y)^2} = \frac{1}{\cos^2(\mathbf s,\mathbf y)}. } \\]

This ratio measures how far the gradient response deviates from being parallel to the displacement.

- Near 1: the observed force response resembles a single effective stiffness direction.
- Much greater than 1: the displacement excites a mixture of modes or strong off-diagonal coupling.

It is not a rigorous condition-number estimate, especially on a changing nonlinear landscape. But it is a cheap, rotation-invariant diagnostic of whether a scalar spectral step is adequate.

I would record it alongside FIRE's velocity–force alignment. It could help decide whether the optimizer should trust a scalar mobility, preserve inertia, or give more weight to a low-rank secant correction.

## 7. What I would actually ask the coding agent to implement

I would avoid another large simultaneous optimizer-development effort. The existing implementation already contains most of the necessary infrastructure.

My suggested experimental order is:

| Priority | Candidate                                                                 | Why test it                                                                         |
| -------- | ------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| 1        | Anderson m=1,2,3, with current safeguards                                 | Determine how much historical rank is actually needed                               |
| 2        | Preconditioned Anderson m=1,2 using the existing particle-secant mobility | Combine diagonal stiffness adaptation with low-rank coupling information            |
| 3        | Persistent adaptive mobility FIRE with constant nominal step              | Test your mass-instead-of-timestep idea directly                                    |
| 4        | FIRE + rank-one stiff-mode suppression                                    | Target coupled modes that per-particle masses cannot isolate                        |
| 5        | Inertial secant hybrid                                                    | Combine useful early FIRE dynamics with progressively reliable Anderson corrections |
| 6        | Tiny learned or symbolic controller                                       | Only after identifying which extra mechanisms actually help                         |

For the first experiment, the existing `DynamicOpt::run_GH()` already supports the appropriate Anderson formulation. Testing m=1, 2 and 3 should mostly require parameter changes, along with checks that the least-squares solver and safeguard behave correctly at those ranks.

For the second, I would use the current secant mass estimates to construct a diagonal initial inverse-Jacobian approximation \\(\eta D\\). It is important to adapt this baseline slowly or freeze it over several accepted steps, since continually changing \\(D\\) makes the secant interpretation less stable.

For the third, use a fixed global nominal scale while adapting the per-particle mobility. Measure how many global stability interventions remain necessary. It would be interesting if a persistent stiffness estimate largely eliminated FIRE's repeated timestep reductions.

### Benchmark design

Run the experiments on both synthetic and physical cases, but keep the first comparison narrow.

For synthetic tests, use rotated anisotropic quadratics, two strongly coupled stiff/soft coordinates embedded in a larger space, and curved valleys with adjustable stiffness and turning radius.

For FireCore, start with the corrected frozen-core H₂O free- and fixed-nuclei cases, then expand to several distorted initial configurations, all-electron eFF, UFF molecules such as xylitol, and RARFF.

Measure actual force-evaluation count, wall time, convergence success, final raw force, reset frequency, and rejected steps. Compare at several force tolerances.

The 72-evaluation BFGS result in the document came from an earlier benchmark at a different convergence threshold, so it should be remeasured under the same final protocol rather than treated as a perfectly matched data point.

I would also make sure the documented remaining energy–force inconsistency in eFF's electron-width derivative is resolved before interpreting failures of line-search-based variants at tight tolerances. The audit identifies this as an outstanding issue.&#x20;

eFF_optimizer_speedup.md



## 8. What I find most promising

I would now distinguish three timescales in our optimization algorithm.

Slow memory — stiffness

Persistent per-particle mobility estimates, updated conservatively. They describe which degrees of freedom are consistently difficult to accelerate.

Intermediate memory — curvature coupling

One or two secant pairs that capture how the gradient responds to recent displacements. Forget or replace them when the local landscape changes.

Fast memory — inertia

Momentum, velocity quenching and short-term steering. These respond immediately to overshooting and changing descent directions.

The important insight is that FIRE and Anderson need not remain separate algorithms.

FIRE's inertia remembers where the optimizer has been moving. Anderson's secant history remembers how the force field responded to that motion. Adaptive masses remember which degrees of freedom have been consistently stiff.

These are complementary forms of memory, operating on different timescales.

I would make the preconditioned one- or two-secant optimizer the immediate next experiment. It is mathematically simple, has a clear relationship to Anderson, requires only a few full-size vectors, and can be implemented using much of what already exists in `DynamicOpt`.

If it captures a substantial fraction of m=8 Anderson's improvement, we can then experiment with mixing it continuously with FIRE's inertia. That seems a more promising route to a genuinely new, compact optimizer than further tuning FIRE's five or six conventional parameters alone.