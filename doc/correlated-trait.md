# The generation of genetically correlated traits
## Correlated traits
Suppose we have a base population with $n$ unrelated ID, and $k$ traits recorded.

Traditionally, it is fairly easy to simulate true breeding values (TBV) of these $k$ traits.  For example, suppose these $k$ traits have zero means and a covariance structure of unit matrix $\mathbf{I}$, then in Julia:

```julia
using Random, Distributions, LinearAlgebra, Statistics
n = 1000
k = 5                  # for example.
tbv = rand(MvNormal(zeros(k), I), n)'
cov(tbv)               # To check if covariance I are realized.
```

Note, $\mathbf{I}$ here can be replaced with an arbitrary positive definite symmetric $\Sigma$.

This is, however, not the case where TBV are simulated with known QTL.

Let 

- $m_i$ be the number of QTL,
- $\mathbf{Q}_i$ be the QTL genotypes of dimension $m_i\times n$,
- $\mathbf{a}_i$ be the QTL effect column vector of length $m_i$,
- $\mathbf{t}_i$ be the TBV column vector of length $n$

of trait $i$.  Then the TBV of trait $i$ is $\mathbf{t}_i=\mathbf{Q}_i'\mathbf{a}$

If we simulate $\mathbf{a}_i$ of each trait independently, then each trait will have an expectation $\mu_i\approx$ mean($\mathbf{t}_i$), and variance $\sigma_i^2\approx$ Var($\mathbf{t}_i$).  Because of the central limit theory, when $m_i$ is large enough, $t_{ij}\stackrel{\mathrm{iid}}{\sim}N(\mu_i, \sigma^2_i)$, for ID $j$ in this population, no matter what distribution of $\mathbf{a}_i$ is.  The off-diagonals of $\Sigma$ of these $k$ traits are unpredictable, as they depend on the distributions of $\mathbf{a}$.  For example, if we independently simulate two traits of the same population with $e_i\sim\Gamma(0.4)$, the two trait will have a correlation of approximately 0.2 on average, and bell-shape distributed.  With different parameters, or distributions, the covariance of the two traits is different.

It is convenient to scale down $\mathbf{a}_{i}$ with $\sigma_i$, such that the covariance matrix $\Sigma$ has unit diagonals, and correlation between traits as off-diagonals.

To construct correlated traits with desired $\rho$ is beyond my ability.  For,

1. Since $\mathbf{a}_i$ of various distribution can all leads to a normal distribution, it is impossible to simulate correlated TBV, then calculate $\mathbf{a}_i$ of desired distribution.
2. A convenient reasoning is to construct a $\mathbf{T}$, such that new $\mathbf{a}_i^* = \mathbf{Ta}_i$ may resulted in a desired $\rho$.  But the linear combination of $\mathbf{a}_i$ has no simple distribution if $\mathbf{a}$ are not from a normal distribution.
3. One brute way is to sample $\mathbf{a}_i$ until $\rho$ is reached.  But this constraint may resulted in a non desired distribution of $\mathbf{a}_i$.  With more than 2 traits, this becomes more difficult.
4. The genetic meaning of such $\mathbf{a}_i$ is not clear, even if they are successfully generated.

My conclusion is it is better to simulate traits with pleiotropic QTL.  The set/number of pleiotropic QTL can be specified, e.g. 10% QTL are shared by two traits.

## Default Laplace and $\Gamma(0.4)$ are very different
![Laplace(0,1) vs $\Gamma(0.4)$](fig/laplace-vs-gamma.pdf)

## One trait
### QTL effect distribution won't affect TBV distribution
Maybe it is because of the central limit theory, but I just show some example here:

In the above figure, I simulated a large population of 10k ID, 24 chromosomes, and 96k SNP of MAF=0.01.

1. Sample 1000/100/10 loci
2. Sample QTL effects
   - $\Gamma(0.4, 1)$
   - `Laplace(0, 1)`
   - $N(0,1)$
   - Uniform(0, 1)
   - Weibull(1, 1)
   - $\chi^2(1)$
3. Calculate TBV
4. Plot nomalized histogram of the TBV

![10 QTL](fig/tbv-10.pdf)

![100 QTL](fig/tbv-100.pdf)

![1000 QTL](fig/tbv-1000.pdf)

To simulate a trait with QTL of extreme values,
it is often to use a $\Gamma(0.4)$ distribution to simulate QTL effect $a_i$.

Let $\mathbf{Q}$ be the QTL genotypes of the 2nd trait of dimension $N_{\mathrm{QTL}}\times N_{\mathrm{ID}}$.  Then TBV$_2=\mathbf{Q'e}$, where $\mathbf{e}$ is the QTL effect vector.
For the first trait $e_i\stackrel{iid}{\sim}\Gamma(0.4)$.

The problem is then reduced to sample $e_i$ such that cov(TBV$_1$, TBV$_2$) equals a given number $\rho$.  We will use $\mathbf{t}_1$ and $\mathbf{t}_2$ below for the 2 TBV.

