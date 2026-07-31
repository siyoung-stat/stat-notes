# Sampling Distributions and Their Approximations

> 송성주·전명식, 『수리통계학』 제5판(자유아카데미, 2020) 3장 학습 노트

## 3.1 Population and Sample

### Definition 3.1.1 — population, sample

- **population**: the totality of objects under study
- **sample**: the part of the population actually observed

### Definition 3.1.2 — random sample

$X_1, \dots, X_n$ is a random sample of size $n$ from population density $f$ if the joint density factors into the product of $f$.

$$f_{X_1,\dots,X_n}(x_1, \dots, x_n) = f(x_1)f(x_2)\cdots f(x_n) = \prod_{i=1}^{n} f(x_i)$$

That is, the $X_i$ are **iid** (independent and identically distributed).

### Example 3.1.1 — Bernoulli population

$f(x) = p^x q^{1-x}$, $x = 0, 1$, $p + q = 1$, $n = 10$:

$$f_{X_1,\dots,X_{10}}(x_1, \dots, x_{10}) = \prod_{i=1}^{10} p^{x_i} q^{1-x_i} = p^{\sum x_i} q^{10 - \sum x_i} \qquad \text{(by independence)}$$

### Definition 3.1.3 — statistic

A **statistic** is a function $T = T(X_1, \dots, X_n)$ of the random sample that contains no unknown parameter. Being a function of random variables, $T$ is itself a random variable.

### Example 3.1.2 — identifying a statistic

A function that depends on an unknown parameter $\theta$ is not a statistic.

- statistic: $\displaystyle \bar{X}_n = \frac{1}{n}\sum_{i=1}^{n} X_i$, $\quad \max\{X_1, \dots, X_n\}$
- not a statistic: $\bar{X}_n - \theta$, $\quad \max\{X_1/\theta, \dots, X_n/\theta\}$

### Definition 3.1.4 — sample moments

- $r$-th sample moment: $\displaystyle m_r' = \frac{1}{n}\sum_{i=1}^{n} X_i^r$
- $r$-th central sample moment: $\displaystyle m_r = \frac{1}{n}\sum_{i=1}^{n} (X_i - \bar{X}_n)^r$

### Definition 3.1.5 — sample mean, sample variance

- sample mean: $\displaystyle \bar{X}_n = \frac{1}{n}\sum_{i=1}^{n} X_i \;\; (= m_1')$
- sample variance: $\displaystyle S_n^2 = \frac{1}{n-1}\sum_{i=1}^{n} (X_i - \bar{X}_n)^2$

Unlike $m_2$, the divisor is $n-1$ rather than $n$; this makes $E(S_n^2) = \sigma^2$, so $S_n^2$ is an unbiased estimator of $\sigma^2$.

### Theorem 3.1.1 — moments of $\bar{X}_n$ and $S_n^2$

For a random sample from a population with mean $\mu$, variance $\sigma^2$, and 4th central moment $\mu_4 = E(X_i - \mu)^4$:

1. $E(\bar{X}_n) = \mu$
2. $\displaystyle \mathrm{Var}(\bar{X}_n) = \frac{\sigma^2}{n}$
3. $E(S_n^2) = \sigma^2$
4. $\displaystyle \mathrm{Var}(S_n^2) = \frac{1}{n}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right), \quad n > 1$

**Proof.**

**(1)**

$$E(\bar{X}_n) = \frac{1}{n}\sum_{i=1}^{n} E(X_i) = \frac{1}{n} n\mu = \mu \qquad \text{(by linearity of } E)$$

**(2)**

$$\mathrm{Var}(\bar{X}_n) = \frac{1}{n^2}\mathrm{Var}\left(\sum_{i=1}^{n} X_i\right) = \frac{1}{n^2}\sum_{i=1}^{n}\mathrm{Var}(X_i) = \frac{\sigma^2}{n} \qquad \text{(by independence)}$$

**(3)** Start from the identity

$$\sum_{i=1}^{n}(X_i - \bar{X}_n)^2 = \sum_{i=1}^{n}(X_i - \mu)^2 - n(\bar{X}_n - \mu)^2$$

Taking expectations on both sides,

$$E\left[\sum_{i=1}^{n}(X_i - \bar{X}_n)^2\right] = n\sigma^2 - n\cdot\frac{\sigma^2}{n} = (n-1)\sigma^2$$

(by $E(X_i - \mu)^2 = \sigma^2$, $E(\bar{X}_n - \mu)^2 = \mathrm{Var}(\bar{X}_n)$)

$$\therefore\; E(S_n^2) = \frac{1}{n-1}(n-1)\sigma^2 = \sigma^2$$

**(4)** Put $Y_i = X_i - \mu$, so $E(Y_i) = 0$, $E(Y_i^2) = \sigma^2$, $E(Y_i^4) = \mu_4$. Since the $Y_i$ are iid, any term in which a distinct index survives with exponent 1 has expectation 0.

$$W := (n-1)S_n^2 = \sum_{i=1}^{n} Y_i^2 - n\bar{Y}_n^2 \quad \text{(by the identity in (3))}, \qquad E(W) = (n-1)\sigma^2$$

$$\mathrm{Var}(S_n^2) = \frac{\mathrm{Var}(W)}{(n-1)^2} = \frac{E(W^2) - [E(W)]^2}{(n-1)^2} = \frac{E(W^2)}{(n-1)^2} - \sigma^4$$

$$W^2 = \left(\sum_{i=1}^{n} Y_i^2\right)^2 - 2n\left(\sum_{i=1}^{n} Y_i^2\right)\bar{Y}_n^2 + n^2\bar{Y}_n^4$$

**(i)** In $E\left[\left(\sum_i Y_i^2\right)^2\right] = \sum_i\sum_j E(Y_i^2 Y_j^2)$, the $n$ terms with $i = j$ give $\mu_4$ and the $n(n-1)$ terms with $i \neq j$ give $\sigma^4$ (by independence):

$$E\left[\left(\sum_{i=1}^{n} Y_i^2\right)^2\right] = n\mu_4 + n(n-1)\sigma^4$$

**(ii)** Substituting $\bar{Y}_n^2 = \frac{1}{n^2}\sum_{j,k} Y_j Y_k$, only the terms with $j = k$ survive (by independence, $E(Y_j) = 0$):

$$E\left[\left(\sum_{i=1}^{n} Y_i^2\right)\bar{Y}_n^2\right] = \frac{1}{n^2}\left[n\mu_4 + n(n-1)\sigma^4\right] = \frac{\mu_4 + (n-1)\sigma^4}{n}$$

**(iii)** In $\bar{Y}_n^4 = \frac{1}{n^4}\sum_{i,j,k,l} Y_i Y_j Y_k Y_l$, only two cases survive: all four indices equal ($n$ terms, $\mu_4$) and two matched pairs ($3n(n-1)$ terms, $\sigma^4$):

$$E(\bar{Y}_n^4) = \frac{1}{n^4}\left[n\mu_4 + 3n(n-1)\sigma^4\right] = \frac{\mu_4 + 3(n-1)\sigma^4}{n^3}$$

Substituting (i), (ii), (iii),

$$E(W^2) = \left[n\mu_4 + n(n-1)\sigma^4\right] - 2n\cdot\frac{\mu_4 + (n-1)\sigma^4}{n} + n^2\cdot\frac{\mu_4 + 3(n-1)\sigma^4}{n^3}$$

Collecting the coefficients of $\mu_4$ and $\sigma^4$,

- coefficient of $\mu_4$: $\displaystyle n - 2 + \frac{1}{n} = \frac{(n-1)^2}{n}$
- coefficient of $\sigma^4$: $\displaystyle n(n-1) - 2(n-1) + \frac{3(n-1)}{n} = \frac{(n-1)(n^2 - 2n + 3)}{n}$

$$E(W^2) = \frac{(n-1)^2}{n}\mu_4 + \frac{(n-1)(n^2 - 2n + 3)}{n}\sigma^4$$

Dividing by $(n-1)^2$, subtracting $\sigma^4$, and putting the $\sigma^4$ term over $n(n-1)$,

$$\frac{n^2 - 2n + 3}{n(n-1)} - 1 = \frac{(n^2 - 2n + 3) - (n^2 - n)}{n(n-1)} = -\frac{n-3}{n(n-1)}$$

$$\therefore\; \mathrm{Var}(S_n^2) = \frac{E(W^2)}{(n-1)^2} - \sigma^4 = \frac{\mu_4}{n} - \frac{n-3}{n(n-1)}\sigma^4 = \frac{1}{n}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right) \qquad \blacksquare$$

---

## 3.2 Distributions Related to the Normal

## 3.2.1 Chi-squared distribution

### Definition 3.2.1 — chi-squared distribution

The gamma distribution with parameters $(k, \theta) = (n/2, 2)$ is called the **chi-squared distribution** with $n$ degrees of freedom, written $X \sim \chi^2(n)$.

$$f_X(x) = \frac{1}{\Gamma(n/2)\,2^{n/2}}\, x^{n/2 - 1} e^{-x/2}, \quad x > 0 \tag{3.1}$$

### Theorem 3.2.1 — mgf, mean, variance of $\chi^2(n)$

1. $M_X(t) = (1 - 2t)^{-n/2}$, $\; t < \tfrac{1}{2}$
2. $E(X) = n$
3. $\mathrm{Var}(X) = 2n$

**Proof.** $\chi^2(n) = \mathrm{GAM}(n/2, 2)$ and the gamma mgf is $(1 - \theta t)^{-k}$, so putting $k = n/2$, $\theta = 2$,

$$M_X(t) = (1 - 2t)^{-n/2} \qquad \text{(by the gamma mgf)}$$

$$E(X) = k\theta = \frac{n}{2}\cdot 2 = n, \qquad \mathrm{Var}(X) = k\theta^2 = \frac{n}{2}\cdot 2^2 = 2n$$

(by the gamma mean and variance) $\blacksquare$

### Theorem 3.2.2 — additivity of the chi-squared distribution

For independent $X_i \sim \chi^2(k_i)$, $i = 1, \dots, n$:

$$Y = \sum_{i=1}^{n} X_i \;\sim\; \chi^2\left(\sum_{i=1}^{n} k_i\right)$$

**Proof.**

$$M_Y(t) = \prod_{i=1}^{n}(1 - 2t)^{-k_i/2} = (1 - 2t)^{-\sum_i k_i/2} \qquad \text{(by the mgf of an independent sum, Theorem 3.2.1)}$$

This is the mgf of $\chi^2\left(\sum_i k_i\right)$, hence the claim (by uniqueness of the mgf). $\blacksquare$

### Theorem 3.2.3 — square of a standard normal

$Z \sim N(0,1) \;\Rightarrow\; Y = Z^2 \sim \chi^2(1)$.

**Proof.**

$$M_{Z^2}(t) = E(e^{tZ^2}) = \int_{-\infty}^{\infty}\frac{1}{\sqrt{2\pi}}e^{tz^2 - z^2/2}\,dz = \int_{-\infty}^{\infty}\frac{1}{\sqrt{2\pi}}e^{-z^2(1-2t)/2}\,dz = (1 - 2t)^{-1/2}, \quad t < \tfrac{1}{2}$$

(by the Gaussian integral). This is the mgf of $\chi^2(1)$, hence the claim (by uniqueness of the mgf). $\blacksquare$

### Theorem 3.2.4 — sum of squared standardized normals

For independent $X_i \sim N(\mu_i, \sigma_i^2)$, $i = 1, \dots, k$:

$$V = \sum_{i=1}^{k}\left(\frac{X_i - \mu_i}{\sigma_i}\right)^2 \;\sim\; \chi^2(k)$$

**Proof.** The standardized variables $(X_i - \mu_i)/\sigma_i \sim N(0,1)$ are independent, each square is $\chi^2(1)$ (by Theorem 3.2.3), and the sum is $\chi^2(k)$ (by Theorem 3.2.2). $\blacksquare$

### Corollary 3.2.1

For a random sample from $N(\mu, \sigma^2)$:

$$\sum_{i=1}^{n}\frac{(X_i - \mu)^2}{\sigma^2} \;\sim\; \chi^2(n)$$

## 3.2.2 F distribution

### Definition 3.2.2 — F distribution

For independent $U \sim \chi^2(n)$ and $V \sim \chi^2(m)$, the distribution of

$$X = \frac{U/n}{V/m} \tag{3.2}$$

is the **F distribution** with degrees of freedom $(n, m)$, written $X \sim F(n, m)$, with density

$$f_X(x) = \frac{\Gamma[(n+m)/2]}{\Gamma(n/2)\Gamma(m/2)}\left(\frac{n}{m}\right)^{n/2}\frac{x^{(n-2)/2}}{[1 + nx/m]^{(n+m)/2}}, \quad x > 0 \tag{3.3}$$

**Derivation of (3.3).** Bivariate change of variables: for a one-to-one map with inverse $u = g(x,y)$, $v = h(x,y)$,

$$f_{X,Y}(x,y) = f_{U,V}\big(g(x,y), h(x,y)\big)\,|J|, \qquad J = \frac{\partial(u,v)}{\partial(x,y)}$$

**(1)** Joint density of $(U, V)$, the product of two densities (3.1) (by independence):

$$f_{U,V}(u,v) = \frac{1}{\Gamma(n/2)\Gamma(m/2)2^{(n+m)/2}}\,u^{n/2-1}v^{m/2-1}e^{-(u+v)/2}, \quad u > 0,\; v > 0$$

**(2)** Attach the auxiliary variable $Y = V$; the inverse map and Jacobian are

$$u = \frac{n}{m}xy, \quad v = y; \qquad J = \begin{vmatrix} \dfrac{n}{m}y & \dfrac{n}{m}x \\[4pt] 0 & 1 \end{vmatrix} = \frac{n}{m}y$$

$$f_{X,Y}(x,y) = \frac{(n/m)^{n/2}}{\Gamma(n/2)\Gamma(m/2)2^{(n+m)/2}}\,x^{n/2-1}y^{(n+m)/2-1}\exp\left[-\frac{y}{2}\left(1 + \frac{nx}{m}\right)\right]$$

(exponent of $y$: $(n/2 - 1) + (m/2 - 1) + 1 = (n+m)/2 - 1$)

**(3)** Marginalize over $y$ with $a = \tfrac{1}{2}(1 + nx/m)$:

$$\int_0^{\infty} y^{(n+m)/2-1}e^{-ay}\,dy = \frac{\Gamma[(n+m)/2]}{a^{(n+m)/2}} = \Gamma[(n+m)/2]\cdot\frac{2^{(n+m)/2}}{(1 + nx/m)^{(n+m)/2}}$$

(by the gamma integral). The factor $2^{(n+m)/2}$ cancels, giving (3.3). $\blacksquare$

### Theorem 3.2.5 — mean and variance of $F(n, m)$

1. $\displaystyle E(X) = \frac{m}{m-2}$, $\; m > 2$
2. $\displaystyle \mathrm{Var}(X) = \frac{2m^2(n + m - 2)}{n(m-2)^2(m-4)}$, $\; m > 4$

**Proof.** $X = \frac{m}{n}\cdot\frac{U}{V}$, where $U$ and $1/V$ are independent and $E(1/V) = 1/(m-2)$ for $m > 2$, so

$$E(X) = \frac{m}{n}E(U)E\!\left(\frac{1}{V}\right) = \frac{m}{n}\cdot n\cdot\frac{1}{m-2} = \frac{m}{m-2} \qquad \text{(by the gamma integral)}$$

For the variance, $X^2 = \dfrac{m^2}{n^2}\cdot\dfrac{U^2}{V^2}$ and

$$E(U^2) = \mathrm{Var}(U) + [E(U)]^2 = 2n + n^2 = n(n+2) \qquad \text{(by Theorem 3.2.1)}$$

$$E\!\left(\frac{1}{V^2}\right) = \int_0^{\infty}\frac{1}{v^2}\cdot\frac{v^{m/2-1}e^{-v/2}}{\Gamma(m/2)2^{m/2}}\,dv = \frac{\Gamma(m/2 - 2)2^{m/2-2}}{\Gamma(m/2)2^{m/2}} = \frac{1}{(m-2)(m-4)}, \quad m > 4$$

(by the gamma integral and $\Gamma(m/2) = (m/2 - 1)(m/2 - 2)\Gamma(m/2 - 2)$)

$$E(X^2) = \frac{m^2}{n^2}E(U^2)E\!\left(\frac{1}{V^2}\right) = \frac{m^2(n+2)}{n(m-2)(m-4)} \qquad \text{(by independence)}$$

$$\mathrm{Var}(X) = E(X^2) - [E(X)]^2 = \frac{m^2(n+2)}{n(m-2)(m-4)} - \frac{m^2}{(m-2)^2} = \frac{2m^2(n + m - 2)}{n(m-2)^2(m-4)}$$

(numerator: $(n+2)(m-2) - n(m-4) = 2(n + m - 2)$) $\blacksquare$

### Theorem 3.2.6 — reciprocal of an F variable, percentile relation

$X \sim F(n, m) \;\Rightarrow\; 1/X \sim F(m, n)$, since $1/X = \dfrac{V/m}{U/n}$. Hence for upper percentiles

$$\frac{1}{F_{\alpha}(n, m)} = F_{1-\alpha}(m, n) \tag{3.4}$$

so the lower part of the F table need not be tabulated.

### Example 3.2.1 — F percentile

$X \sim F(5, 10)$ and $P(X \leq a) = 0.01$:

$$a = F_{0.99}(5,10) = \frac{1}{F_{0.01}(10,5)} = \frac{1}{10.1} = 0.099 \qquad \text{(by (3.4) and the F table)}$$

## 3.2.3 t distribution

### Definition 3.2.3 — t distribution

For $Z \sim N(0,1)$ and independent $U \sim \chi^2(k)$, the distribution of

$$X = \frac{Z}{\sqrt{U/k}} \tag{3.5}$$

is the **t distribution** with $k$ degrees of freedom, written $X \sim t(k)$, with density

$$f_X(x) = \frac{\Gamma[(k+1)/2]}{\Gamma(k/2)\sqrt{k\pi}}\cdot\frac{1}{(1 + x^2/k)^{(k+1)/2}}, \quad -\infty < x < \infty \tag{3.6}$$

**Derivation of (3.6).** The same bivariate change of variables, with auxiliary variable $W = U$.

**(1)** Joint density (by independence):

$$f_{Z,U}(z,u) = \frac{1}{\sqrt{2\pi}}e^{-z^2/2}\cdot\frac{1}{\Gamma(k/2)2^{k/2}}u^{k/2-1}e^{-u/2}, \quad u > 0$$

**(2)** Inverse map and Jacobian:

$$z = x\sqrt{\frac{w}{k}}, \quad u = w; \qquad J = \begin{vmatrix} \dfrac{\sqrt{w}}{\sqrt{k}} & \dfrac{x}{2\sqrt{kw}} \\[4pt] 0 & 1 \end{vmatrix} = \sqrt{\frac{w}{k}}$$

$$f_{X,W}(x,w) = \frac{1}{\sqrt{k\pi}\,\Gamma(k/2)2^{(k+1)/2}}\,w^{(k+1)/2-1}\exp\left[-\frac{w}{2}\left(1 + \frac{x^2}{k}\right)\right]$$

(exponent of $w$: $w^{k/2-1}\cdot w^{1/2} = w^{(k+1)/2-1}$; constants: $\sqrt{2\pi}\cdot 2^{k/2}\cdot\sqrt{k} = 2^{(k+1)/2}\sqrt{k\pi}$)

**(3)** Marginalize over $w$ with $a = \tfrac{1}{2}(1 + x^2/k)$:

$$\int_0^{\infty} w^{(k+1)/2-1}e^{-aw}\,dw = \Gamma[(k+1)/2]\cdot\frac{2^{(k+1)/2}}{(1 + x^2/k)^{(k+1)/2}} \qquad \text{(by the gamma integral)}$$

The factor $2^{(k+1)/2}$ cancels, giving (3.6). $\blacksquare$

### Example 3.2.2 — symmetry of the t distribution

$X \sim t(10)$ with $P(X \leq 1.812) = 0.95$ from the table; the density is symmetric about $x = 0$, so

$$P(X \leq -1.812) = 1 - P(X \leq 1.812) = 0.05$$

### Theorem 3.2.7 — mean and variance of $t(n)$

1. $E(X) = 0$, $\; n > 1$
2. $\displaystyle \mathrm{Var}(X) = \frac{n}{n-2}$, $\; n > 2$

**Proof.** $Z$ is symmetric about 0 and independent of $U$, so $X$ is symmetric about 0; for $n > 1$, $E|X| < \infty$ and hence $E(X) = 0$. Then $\mathrm{Var}(X) = E(X^2)$ and

$$X^2 = \frac{Z^2}{U/n} = nZ^2\cdot\frac{1}{U}, \qquad E(X^2) = nE(Z^2)E\!\left(\frac{1}{U}\right) \quad \text{(by independence)}$$

$$E(Z^2) = \mathrm{Var}(Z) + [E(Z)]^2 = 1$$

$$E\!\left(\frac{1}{U}\right) = \frac{1}{\Gamma(n/2)2^{n/2}}\int_0^{\infty} u^{(n/2-1)-1}e^{-u/2}\,du = \frac{\Gamma(n/2 - 1)2^{n/2-1}}{\Gamma(n/2)2^{n/2}} = \frac{1}{n-2}, \quad n > 2$$

(by the gamma integral $\int_0^{\infty} u^{a-1}e^{-u/2}\,du = \Gamma(a)2^a$ and $\Gamma(n/2) = (n/2 - 1)\Gamma(n/2 - 1)$)

$$\therefore\; \mathrm{Var}(X) = E(X^2) = n\cdot 1\cdot\frac{1}{n-2} = \frac{n}{n-2} \qquad \blacksquare$$

### Corollary 3.2.2 — square of a t variable

$X \sim t(n) \;\Rightarrow\; X^2 \sim F(1, n)$.

**Proof.**

$$X^2 = \frac{Z^2}{U/n} = \frac{Z^2/1}{U/n}$$

with $Z^2 \sim \chi^2(1)$ (by Theorem 3.2.3) independent of $U \sim \chi^2(n)$, which is exactly the definition (3.2) of $F(1, n)$. $\blacksquare$

### Notation 3.2.1 — upper percentiles

The point with upper-tail probability $\alpha$ is written as follows.

- $z_{\alpha}$: $P(Z > z_{\alpha}) = \alpha$ for $Z \sim N(0,1)$
- $t_{\alpha}(n)$: $P(T > t_{\alpha}(n)) = \alpha$ for $T \sim t(n)$
- $\chi_{\alpha}^2(n)$: $P(X > \chi_{\alpha}^2(n)) = \alpha$ for $X \sim \chi^2(n)$
- $F_{\alpha}(n, m)$: $P(F > F_{\alpha}(n,m)) = \alpha$ for $F \sim F(n,m)$

$N(0,1)$ and $t(n)$ are symmetric about 0, hence

$$z_{1-\alpha} = -z_{\alpha}, \qquad t_{1-\alpha}(n) = -t_{\alpha}(n)$$

---

## 3.3 Sampling from a Normal Distribution

### Theorem 3.3.1 — sum of independent normals

For independent $X_i \sim N(\mu_i, \sigma_i^2)$:

$$\sum_{i=1}^{n} X_i \;\sim\; N\left(\sum_{i=1}^{n}\mu_i,\; \sum_{i=1}^{n}\sigma_i^2\right)$$

**Proof.** The normal mgf is $\exp(\mu_i t + \sigma_i^2 t^2/2)$, so

$$M(t) = \prod_{i=1}^{n}\exp\left(\mu_i t + \frac{\sigma_i^2 t^2}{2}\right) = \exp\left[t\sum_{i=1}^{n}\mu_i + \frac{t^2}{2}\sum_{i=1}^{n}\sigma_i^2\right]$$

(by the mgf of an independent sum). This is the mgf of $N\left(\sum\mu_i, \sum\sigma_i^2\right)$, hence the claim (by uniqueness of the mgf). $\blacksquare$

### Corollary 3.3.1 — distribution of the sample mean

For a random sample of size $n$ from $N(\mu, \sigma^2)$:

$$\bar{X}_n \;\sim\; N\left(\mu, \frac{\sigma^2}{n}\right)$$

(by Theorem 3.3.1 and the linear transformation property of the normal)

### Example 3.3.1 — production time

$X \sim N(6, 2^2)$ independent, $n = 10$; then $\sum X_i \sim N(60, 40)$ (by Theorem 3.3.1) and

$$P\left(\sum_{i=1}^{10} X_i \geq 70\right) = P\left[\frac{\sum X_i - 60}{\sqrt{40}} \geq \frac{70 - 60}{\sqrt{40}}\right] = 1 - \Phi(1.58) = 0.057$$

### Theorem 3.3.2 — independence of $\bar{X}_n$ and $S_n^2$, and the distribution of $S_n^2$

For a random sample of size $n$ from $N(\mu, \sigma^2)$:

1. $\bar{X}_n$ and $S_n^2$ are independent
2. $\displaystyle \frac{(n-1)S_n^2}{\sigma^2} \sim \chi^2(n-1)$

Note that the deviations $X_i - \bar{X}_n$ are **not** independent, because of the constraint $\sum_i (X_i - \bar{X}_n) = 0$; Theorem 3.2.4 therefore cannot be applied to $S_n^2$ directly.

**Proof of (1), by an orthogonal transformation.** An orthogonal matrix satisfies $Q^{\top}Q = QQ^{\top} = I$ and preserves squared length:

$$\|Qx\|^2 = x^{\top}Q^{\top}Qx = x^{\top}x = \|x\|^2$$

For an iid sample the covariance matrix is $\mathrm{Cov}(X) = \sigma^2 I$, so for $Y = QX$

$$\mathrm{Cov}(Y) = Q\mathrm{Cov}(X)Q^{\top} = Q(\sigma^2 I)Q^{\top} = \sigma^2 QQ^{\top} = \sigma^2 I$$

$Y$ is multivariate normal (a linear transformation of a multivariate normal), and with $\Sigma = \sigma^2 I$ the quadratic form in the exponent separates componentwise, so the joint density factors into its marginals: $Y_1, \dots, Y_n$ are mutually independent, and functions of disjoint index sets are independent.

Take $Q$ whose first row is $(1/\sqrt{n}, \dots, 1/\sqrt{n})$ (a **Helmert matrix**). Then

$$Y_1 = \frac{1}{\sqrt{n}}\sum_{i=1}^{n} X_i = \sqrt{n}\,\bar{X}_n$$

$$\sum_{i=2}^{n} Y_i^2 = \sum_{i=1}^{n} X_i^2 - Y_1^2 = \sum_{i=1}^{n} X_i^2 - n\bar{X}_n^2 = \sum_{i=1}^{n}(X_i - \bar{X}_n)^2 = (n-1)S_n^2$$

(by length preservation)

Thus $\bar{X}_n = Y_1/\sqrt{n}$ is a function of $Y_1$ alone and $S_n^2 = \frac{1}{n-1}\sum_{i=2}^{n} Y_i^2$ is a function of $Y_2, \dots, Y_n$ alone; the index sets $\{1\}$ and $\{2, \dots, n\}$ are disjoint, so the two are independent. $\blacksquare$

**Proof of (2), by mgf.** From the identity $\sum(X_i - \mu)^2 = \sum(X_i - \bar{X}_n)^2 + n(\bar{X}_n - \mu)^2$, dividing by $\sigma^2$,

$$V_1 := \sum_{i=1}^{n}\frac{(X_i - \mu)^2}{\sigma^2} = \frac{(n-1)S_n^2}{\sigma^2} + \frac{n(\bar{X}_n - \mu)^2}{\sigma^2}$$

Write $V_2 = \dfrac{(n-1)S_n^2}{\sigma^2}$ and $V_3 = \dfrac{n(\bar{X}_n - \mu)^2}{\sigma^2}$, so that $V_1 = V_2 + V_3$.

- $V_1 \sim \chi^2(n)$ (by Theorem 3.2.4)
- $V_3 = \left((\bar{X}_n - \mu)/(\sigma/\sqrt{n})\right)^2 \sim \chi^2(1)$ (by Corollary 3.3.1 and Theorem 3.2.3)
- $V_2$ is a function of $S_n^2$ and $V_3$ of $\bar{X}_n$, so $V_2 \perp V_3$ (by part (1))

$$M_{V_1}(t) = M_{V_2}(t)M_{V_3}(t) \qquad \text{(by the mgf of an independent sum)}$$

$$M_{V_2}(t) = \frac{M_{V_1}(t)}{M_{V_3}(t)} = \frac{(1 - 2t)^{-n/2}}{(1 - 2t)^{-1/2}} = (1 - 2t)^{-(n-1)/2}, \quad t < \tfrac{1}{2} \qquad \text{(by Theorem 3.2.1)}$$

This is the mgf of $\chi^2(n-1)$, hence $V_2 \sim \chi^2(n-1)$ (by uniqueness of the mgf). $\blacksquare$

Only the normal distribution keeps its components independent under rotation, so the independence of $\bar{X}_n$ and $S_n^2$ is a **characterization** of the normal distribution.

### Example 3.3.2 — probability for the sample variance

Continuing Example 3.3.1 ($n = 10$, $\sigma^2 = 4$), $9S_{10}^2/4 \sim \chi^2(9)$ (by Theorem 3.3.2(2)) and

$$P(S_{10}^2 > 5) = P\left[\frac{9S_{10}^2}{4} > \frac{45}{4}\right] = P[\chi^2(9) > 11.25] \approx 0.259$$

### Theorem 3.3.3 — ratio of two sample variances

Independent random samples of size $n$ from $N(\mu_X, \sigma_X^2)$ and of size $m$ from $N(\mu_Y, \sigma_Y^2)$:

$$F = \frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2} \;\sim\; F(n-1, m-1)$$

**Proof.** $(n-1)S_X^2/\sigma_X^2 \sim \chi^2(n-1)$ and $(m-1)S_Y^2/\sigma_Y^2 \sim \chi^2(m-1)$ (by Theorem 3.3.2(2)), and they are independent because the two samples are. Dividing each by its degrees of freedom,

$$F = \frac{[(n-1)S_X^2/\sigma_X^2]/(n-1)}{[(m-1)S_Y^2/\sigma_Y^2]/(m-1)} = \frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2} \qquad \text{(by Definition 3.2.2)} \quad \blacksquare$$

In particular, if $\sigma_X^2 = \sigma_Y^2$ then $S_X^2/S_Y^2 \sim F(n-1, m-1)$, which is the basis of the test for equality of two variances.

### Theorem 3.3.4 — studentized t statistic

For a random sample from $N(\mu, \sigma^2)$, with $S_n = \sqrt{S_n^2}$:

$$T = \frac{\sqrt{n}(\bar{X}_n - \mu)}{\sqrt{\sum_{i=1}^{n}(X_i - \bar{X}_n)^2/(n-1)}} = \frac{\bar{X}_n - \mu}{S_n/\sqrt{n}} \;\sim\; t(n-1)$$

**Proof.**

$$Z = \frac{\sqrt{n}(\bar{X}_n - \mu)}{\sigma} \sim N(0,1), \qquad U = \frac{(n-1)S_n^2}{\sigma^2} \sim \chi^2(n-1)$$

and $Z, U$ are independent (by Theorem 3.3.2), so

$$T = \frac{Z}{\sqrt{U/(n-1)}} = \frac{\sqrt{n}(\bar{X}_n - \mu)/\sigma}{\sqrt{[(n-1)S_n^2/\sigma^2]/(n-1)}} = \frac{\sqrt{n}(\bar{X}_n - \mu)}{S_n} \;\sim\; t(n-1)$$

(by Definition 3.2.3) $\blacksquare$

**Studentization** replaces the nuisance parameter $\sigma^2$ by $S_n^2$, thereby eliminating it. As $n \to \infty$ the denominator converges to $\sigma$, so $t(n-1) \to N(0,1)$ (by Slutsky's theorem; see Example 3.4.5).

---

## 3.4 Law of Large Numbers and Central Limit Theorem

## 3.4.1 Law of large numbers

### Definition 3.4.1 — convergence in probability

$$\lim_{n\to\infty} P(|X_n - X| \geq \epsilon) = 0 \;\text{ for every } \epsilon > 0 \quad \Longleftrightarrow \quad X_n \xrightarrow{p} X$$

equivalently $\lim_{n\to\infty} P(|X_n - X| < \epsilon) = 1$.

### Theorem 3.4.1 — law of large numbers

For a random sample from a population with mean $\mu < \infty$:

$$\bar{X}_n \xrightarrow{p} \mu$$

**Proof.** Assume $\sigma^2 < \infty$. For any $\epsilon > 0$,

$$P[|\bar{X}_n - \mu| < \epsilon] \geq 1 - \frac{E(\bar{X}_n - \mu)^2}{\epsilon^2} = 1 - \frac{\sigma^2/n}{\epsilon^2} \longrightarrow 1 \quad (n \to \infty)$$

(by Chebyshev's inequality and $E(\bar{X}_n - \mu)^2 = \mathrm{Var}(\bar{X}_n) = \sigma^2/n$, Theorem 3.1.1). $\blacksquare$

### Theorem 3.4.2 — convergence in probability of the sample variance

If $E(X_1^4) < \infty$, then

$$S_n^2 = \frac{1}{n-1}\sum_{i=1}^{n}(X_i - \bar{X}_n)^2 \xrightarrow{p} \sigma^2$$

**Proof.** From Theorem 3.1.1, $E(S_n^2) = \sigma^2$ and $\mathrm{Var}(S_n^2) = \frac{1}{n}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right)$, and $E(X_1^4) < \infty$ means $\mu_4 < \infty$. For any $\epsilon > 0$,

$$P[|S_n^2 - \sigma^2| \geq \epsilon] \leq \frac{\mathrm{Var}(S_n^2)}{\epsilon^2} = \frac{1}{n\epsilon^2}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right) \longrightarrow 0 \qquad \text{(by Chebyshev's inequality)}$$

since the bracket converges to $\mu_4 - \sigma^4$ while the factor $1/n$ vanishes. $\blacksquare$

### Example 3.4.1 — convergence of the sample proportion

For independent Bernoulli($p$) variables, $E(X_1) = p$, so

$$\hat{p}_n = \frac{1}{n}\sum_{i=1}^{n} X_i \xrightarrow{p} p \qquad \text{(by Theorem 3.4.1)}$$

## 3.4.2 Central limit theorem

### Definition 3.4.2 — convergence in distribution

$$\lim_{n\to\infty} F_{X_n}(x) = F_X(x) \;\text{ at every continuity point of } F_X \quad \Longleftrightarrow \quad X_n \xrightarrow{d} X$$

What converges is the **cdf**, not the random variable, and only at continuity points of $F_X$.

### Theorem 3.4.3 — continuity theorem

If $\lim_{n\to\infty} M_n(t) = M(t)$ on an open interval $-h < t < h$ and $M(t)$ is the mgf of a distribution with cdf $F$, then $\lim_{n\to\infty} F_n(x) = F(x)$ at every continuity point of $F$.

For discrete variables the pgf may be used instead: if $\lim_n G_n(s) = G(s)$ for $0 \leq s < 1$, then $\lim_n P(X_n = k) = a_k$ where $G(s) = \sum_k a_k s^k$.

### Theorem 3.4.4 — central limit theorem

For a random sample from a population with mean $\mu$ and variance $\sigma^2 < \infty$:

$$Z_n = \frac{\sum_{i=1}^{n} X_i - E\left(\sum_{i=1}^{n} X_i\right)}{\sqrt{\mathrm{Var}\left(\sum_{i=1}^{n} X_i\right)}} = \frac{\sum_{i=1}^{n}(X_i - \mu)}{\sqrt{n}\,\sigma} \xrightarrow{d} N(0,1)$$

### Lemma 3.4.1 — Lagrange form of the first-order remainder

Let $f$ be twice differentiable on an open interval containing $a$, and define the remainder as everything the first-order approximation misses:

$$R(x) := f(x) - f(a) - f'(a)(x - a) \quad \Longleftrightarrow \quad f(x) = f(a) + f'(a)(x-a) + R(x)$$

Then for some $\xi$ between $a$ and $x$,

$$R(x) = \frac{f''(\xi)}{2}(x - a)^2$$

**Proof.** Fix $x \neq a$ and define $M$ by $R(x) = \dfrac{(x-a)^2}{2}M$. Put

$$g(u) = f(u) - f(a) - f'(a)(u - a) - \frac{(u-a)^2}{2}M$$

so that $g(a) = 0$, $g'(a) = 0$, and $g(x) = 0$ (the last by the definition of $M$).

- $g(a) = g(x) = 0 \Rightarrow g'(\xi_1) = 0$ for some $\xi_1 \in (a, x)$ (by Rolle's theorem)
- $g'(a) = g'(\xi_1) = 0 \Rightarrow g''(\xi_2) = 0$ for some $\xi_2 \in (a, \xi_1)$ (by Rolle's theorem)

Differentiating twice kills the first-order part and leaves $M$:

$$g''(u) = f''(u) - M, \qquad g''(\xi_2) = f''(\xi_2) - M = 0 \iff M = f''(\xi_2)$$

Writing $\xi = \xi_2$ gives the claim. $\blacksquare$

The point used below is that $\xi$ is trapped between $a$ and $x$, so $x \to a$ forces $\xi \to a$.

**Proof of Theorem 3.4.4 (assuming the mgf exists).**

**(1)** The mgf of the centered variable is common to all $i$:

$$m(t) = E\left[e^{t(X_i - \mu)}\right], \quad m(0) = 1, \quad m'(0) = E(X_i - \mu) = 0, \quad m''(0) = E(X_i - \mu)^2 = \sigma^2$$

**(2)** Expanding about 0 and splitting off the leading term,

$$m(t) = m(0) + m'(0)t + \frac{m''(\xi)t^2}{2} = 1 + \frac{\sigma^2 t^2}{2} + \frac{[m''(\xi) - \sigma^2]t^2}{2}, \qquad 0 < \xi < t$$

(by Taylor series with the Lagrange remainder, Lemma 3.4.1)

**(3)** Substituting $t/(\sqrt{n}\sigma)$,

$$M_{Z_n}(t) = \prod_{i=1}^{n} M_{X_i - \mu}\left(\frac{t}{\sqrt{n}\sigma}\right) = \left[m\left(\frac{t}{\sqrt{n}\sigma}\right)\right]^n = \left[1 + \frac{t^2}{2n} + \frac{[m''(\xi) - \sigma^2]t^2}{2n\sigma^2}\right]^n, \qquad 0 < \xi < \frac{t}{\sqrt{n}\sigma}$$

(by the mgf of an independent sum)

**(4)** For fixed $t$, $t/(\sqrt{n}\sigma) \to 0$ and $\xi$ is squeezed to 0; $m''$ is continuous at 0, so $m''(\xi) - \sigma^2 \to m''(0) - \sigma^2 = 0$. The bracket is $1 + (t^2/2 + o(1))/n$, hence

$$\lim_{n\to\infty} M_{Z_n}(t) = \exp(t^2/2) \qquad \text{(by the exponential limit } (1 + c/n)^n \to e^c)$$

**(5)** This is the mgf of $N(0,1)$, so $Z_n \xrightarrow{d} N(0,1)$ (by Theorem 3.4.3 and uniqueness of the mgf). $\blacksquare$

The CLT requires only a finite mean and variance, whatever the shape of the population; it fails for e.g. the Cauchy distribution, which has no mean. Equivalently

$$Z_n = \frac{\bar{X}_n - E(\bar{X}_n)}{\sqrt{\mathrm{Var}(\bar{X}_n)}} = \frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}} \tag{3.8}$$

so the CLT approximates the distribution of $\bar{X}_n$ or of $\sum X_i$ by a normal distribution.

### Example 3.4.2 — normal approximation, continuous case

For a random sample from $U(0,1)$, $E(X_i) = \tfrac{1}{2}$ and $\mathrm{Var}(X_i) = \tfrac{1}{12}$, so $\sum X_i \approx N(n/2,\; n/12)$ and

$$P\left[a \leq \sum_{i=1}^{n} X_i \leq b\right] \approx \Phi\left[\frac{b - n/2}{\sqrt{n/12}}\right] - \Phi\left[\frac{a - n/2}{\sqrt{n/12}}\right] \qquad \text{(by Theorem 3.4.4)}$$

### Example 3.4.3 — normal approximation, binomial case

$X \sim B(n, p)$ has the same distribution as $\sum_{i=1}^{n} X_i$ for independent Bernoulli($p$) variables, with $E(X_i) = p$ and $\mathrm{Var}(X_i) = pq$, $q = 1 - p$. Hence $\left(\sum X_i - np\right)/\sqrt{npq}$ is approximately $N(0,1)$ and

$$P[a \leq X \leq b] \approx \Phi\left[\frac{b - np}{\sqrt{npq}}\right] - \Phi\left[\frac{a - np}{\sqrt{npq}}\right] \qquad \text{(by Theorem 3.4.4)}$$

### Definition 3.4.3 — continuity correction

$P(X = k)$ is the area of a bar of width 1 centered at $k$, whereas the normal approximation integrates a curve; taking the integration limits at the integers $a, b$ (the bar centers) cuts off the outer halves of the two end bars. Extending each limit by $\tfrac{1}{2}$ covers the whole bars:

$$P\left(a - \tfrac{1}{2} \leq X \leq b + \tfrac{1}{2}\right)$$

$$P[a \leq X \leq b] \approx \Phi\left[\frac{b + 1/2 - np}{\sqrt{npq}}\right] - \Phi\left[\frac{a - 1/2 - np}{\sqrt{npq}}\right]$$

### Example 3.4.4 — effect of the continuity correction

$n = 50$, $p = 0.75$, $a = 36$, $b = 37$, so $np = 37.5$ and $npq = 9.375$.

- without correction: $\Phi\left[\dfrac{37 - 37.5}{\sqrt{9.375}}\right] - \Phi\left[\dfrac{36 - 37.5}{\sqrt{9.375}}\right] = 0.1230$
- with correction: $\Phi\left[\dfrac{37.5 - 37.5}{\sqrt{9.375}}\right] - \Phi\left[\dfrac{35.5 - 37.5}{\sqrt{9.375}}\right] = 0.2432$
- exact binomial: $P(36 \leq X \leq 37) = 0.2371$

## 3.4.3 Slutsky's theorem and the delta method

### Theorem 3.4.5 — Slutsky's theorem

If $X_n \xrightarrow{p} c$ for a constant $c$ and $Y_n \xrightarrow{d} Z$, then

1. $Y_n + X_n \xrightarrow{d} Z + c$
2. $X_n Y_n \xrightarrow{d} cZ$

**Proof of (1) (outline).** For any $\epsilon > 0$, splitting on the event $|X_n - c| < \epsilon$,

$$F_{Y_n + X_n}(z) = P(Y_n + X_n \leq z,\; |X_n - c| < \epsilon) + P(Y_n + X_n \leq z,\; |X_n - c| \geq \epsilon)$$

$$F_{Y_n + X_n}(z) \leq P(Y_n \leq z - c + \epsilon) + P(|X_n - c| \geq \epsilon)$$

and the same splitting applied to $P(Y_n + c + \epsilon \leq z)$ gives the lower bound, so

$$P(Y_n \leq z - c - \epsilon) - P(|X_n - c| \geq \epsilon) \;\leq\; F_{Y_n + X_n}(z) \;\leq\; P(Y_n \leq z - c + \epsilon) + P(|X_n - c| \geq \epsilon)$$

This holds for every $\epsilon > 0$, and both bounds converge to $P(Z \leq z - c) = P(Z + c \leq z)$ wherever $Z$ is continuous at $z - c$ (by $Y_n \xrightarrow{d} Z$ and $X_n \xrightarrow{p} c$). $\blacksquare$

### Example 3.4.5 — limiting distribution of the studentized statistic

For a random sample from a population with mean $\mu$ and variance $\sigma^2 < \infty$ (no normality assumed): $S_n^2 \xrightarrow{p} \sigma^2$ (by Theorem 3.4.2), hence $S_n \xrightarrow{p} \sigma$ (since $\sqrt{\cdot}$ is continuous), and $(\bar{X}_n - \mu)/(\sigma/\sqrt{n}) \xrightarrow{d} N(0,1)$ (by Theorem 3.4.4). Therefore

$$\frac{\bar{X}_n - \mu}{S_n/\sqrt{n}} = \frac{\sigma}{S_n}\cdot\frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}} \xrightarrow{d} N(0,1) \qquad \text{(by Theorem 3.4.5(2))}$$

Combined with Theorem 3.3.4 this confirms that $t(n-1) \to N(0,1)$ as the degrees of freedom grow.

### Theorem 3.4.6 — delta method

If $\sqrt{n}(\bar{X}_n - \theta) \xrightarrow{d} N(0, \sigma^2)$ and $g'$ is continuous and nonzero at $\theta$, then

$$\sqrt{n}\left(g(\bar{X}_n) - g(\theta)\right) \xrightarrow{d} N\!\left(0,\; \sigma^2[g'(\theta)]^2\right)$$

**Proof.** For some $\tilde{\theta}$ between $\bar{X}_n$ and $\theta$,

$$g'(\tilde{\theta}) = \frac{g(\bar{X}_n) - g(\theta)}{\bar{X}_n - \theta} \quad \Longleftrightarrow \quad g(\bar{X}_n) = g'(\tilde{\theta})(\bar{X}_n - \theta) + g(\theta) \qquad \text{(by the mean value theorem)}$$

$$\bar{X}_n - \theta = \frac{1}{\sqrt{n}}\cdot\sqrt{n}(\bar{X}_n - \theta) \xrightarrow{d} 0\cdot N(0,\sigma^2) = 0 \qquad \text{(by Theorem 3.4.5(2))}$$

Convergence in distribution to a constant is convergence in probability, so $\bar{X}_n \xrightarrow{p} \theta$; $\tilde{\theta}$ is trapped between $\bar{X}_n$ and $\theta$, so $\tilde{\theta} \xrightarrow{p} \theta$ (by the squeeze theorem), and $g'(\tilde{\theta}) \xrightarrow{p} g'(\theta)$ since $g'$ is continuous. Hence

$$\sqrt{n}\left(g(\bar{X}_n) - g(\theta)\right) = g'(\tilde{\theta})\cdot\sqrt{n}(\bar{X}_n - \theta) \xrightarrow{d} g'(\theta)N(0,\sigma^2) = N\!\left(0, \sigma^2[g'(\theta)]^2\right)$$

(by Theorem 3.4.5(2)) $\blacksquare$

### Example 3.4.6 — variance stabilizing transformation

For a random sample from Poisson($\lambda$), $E(X_i) = \mathrm{Var}(X_i) = \lambda$, so

$$Z_n = \sqrt{n}(\bar{X}_n - \lambda) \xrightarrow{d} N(0, \lambda) \qquad \text{(by Theorem 3.4.4)}$$

With $g(x) = \sqrt{x}$, $g'(\lambda) = \dfrac{1}{2\sqrt{\lambda}}$ and $[g'(\lambda)]^2 = \dfrac{1}{4\lambda}$, so

$$\sqrt{n}\left(\sqrt{\bar{X}_n} - \sqrt{\lambda}\right) \xrightarrow{d} N\!\left(0, \tfrac{1}{4}\right) \qquad \text{(by Theorem 3.4.6)}$$

The limiting variance no longer depends on $\lambda$: this is a **variance stabilizing transformation**.

### Theorem 3.4.7 — multivariate delta method: approximate mean and variance

For a random vector $X = (X_1, \dots, X_k)$ with $E(X) = \theta = (\theta_1, \dots, \theta_k)$,

$$g(x) \approx g(\theta) + \sum_{i=1}^{k} g_i'(\theta)(x_i - \theta_i), \qquad g_i'(\theta) = \left.\frac{\partial}{\partial x_i}g(x)\right|_{x = \theta} \qquad \text{(by Taylor series)}$$

$$E\,g(X) \approx g(\theta) + \sum_{i=1}^{k} g_i'(\theta)E(X_i - \theta_i) = g(\theta) \qquad \text{(by } E(X_i - \theta_i) = 0) \tag{3.9}$$

$$\mathrm{Var} g(X) \approx \sum_{i=1}^{k}[g_i'(\theta)]^2\mathrm{Var}(X_i) + 2\sum_{i > j} g_i'(\theta)g_j'(\theta)\mathrm{Cov}(X_i, X_j) \tag{3.10}$$

(3.10) follows by expanding $\mathrm{Var} g(X) \approx E\left[\left(\sum_i g_i'(\theta)(X_i - \theta_i)\right)^2\right]$.

### Example 3.4.7 — ratio-type estimator

Estimating $g(\mu_X, \mu_Y) = \dfrac{\mu_X}{\mu_Y}$ by $\dfrac{\bar{X}_n}{\bar{Y}_m}$, the partial derivatives are

$$\frac{\partial}{\partial \mu_X}g = \frac{1}{\mu_Y}, \qquad \frac{\partial}{\partial \mu_Y}g = -\frac{\mu_X}{\mu_Y^2}$$

$$E\left(\frac{\bar{X}_n}{\bar{Y}_m}\right) \approx \frac{\mu_X}{\mu_Y} \qquad \text{(by (3.9))}$$

$$\mathrm{Var}\left(\frac{\bar{X}_n}{\bar{Y}_m}\right) \approx \frac{1}{\mu_Y^2}\mathrm{Var}(\bar{X}_n) + \frac{\mu_X^2}{\mu_Y^4}\mathrm{Var}(\bar{Y}_m) - 2\frac{\mu_X}{\mu_Y^3}\mathrm{Cov}(\bar{X}_n, \bar{Y}_m) \qquad \text{(by (3.10))}$$

$$= \left(\frac{\mu_X}{\mu_Y}\right)^2\left(\frac{\mathrm{Var}(\bar{X}_n)}{\mu_X^2} + \frac{\mathrm{Var}(\bar{Y}_m)}{\mu_Y^2} - 2\frac{\mathrm{Cov}(\bar{X}_n, \bar{Y}_m)}{\mu_X \mu_Y}\right)$$

---

## 3.5 Order Statistics

### Definition 3.5.1 — order statistics

The random sample $X_1, \dots, X_n$ arranged in increasing order:

$$X_{(1)} \leq X_{(2)} \leq \cdots \leq X_{(n)} \tag{3.11}$$

### Theorem 3.5.1 — joint density of the order statistics

For continuous $f$,

$$g(x_{(1)}, \dots, x_{(n)}) = \begin{cases} n!\, f(x_{(1)})f(x_{(2)})\cdots f(x_{(n)}), & x_{(1)} < x_{(2)} < \cdots < x_{(n)} \\[4pt] 0, & \text{otherwise} \end{cases} \tag{3.12}$$

**Proof.** Sorting is not one-to-one: $n!$ different samples map to the same ordered vector, e.g. all $3! = 6$ permutations of $(5, 1, 3)$ map to $(1, 3, 5)$, so the change-of-variable formula $f_Y(y) = f_X(g^{-1}(y))|J|$ cannot be applied directly.

Partition the sample space into the $n!$ disjoint regions on which the ordering is fixed, e.g. $A_1 = \{x_1 < x_2 < \cdots < x_n\}$; on each region sorting is one-to-one and merely permutes coordinates, so $|J_i| = 1$ and the density there is $f(x_{(1)})\cdots f(x_{(n)})$.

For an iid sample the density is invariant under permutation of its arguments, so all $n!$ regions contribute equally; summing them multiplies one region's density by $n!$, which gives (3.12). The boundaries $x_i = x_j$ have probability 0 in the continuous case, so strict and weak inequalities give the same density. $\blacksquare$

### Example 3.5.1 — joint density for $f(x) = 3x^2$

$X_1, X_2, X_3$ a random sample from $f(x) = 3x^2$, $0 < x < 1$; for $0 < x_{(1)} < x_{(2)} < x_{(3)} < 1$,

$$g(x_{(1)}, x_{(2)}, x_{(3)}) = 3!\,(3x_{(1)}^2)(3x_{(2)}^2)(3x_{(3)}^2) = 162\,x_{(1)}^2 x_{(2)}^2 x_{(3)}^2$$

and 0 otherwise.

### Theorem 3.5.2 — density of the $k$-th order statistic

For a population with density $f$ and cdf $F$, with $f(x) > 0$ on $a < x < b$:

$$f_{X_{(k)}}(x) = \frac{n!}{(k-1)!\,(n-k)!}[F(x)]^{k-1}[1 - F(x)]^{n-k}f(x), \quad a < x < b$$

and 0 otherwise.

**Proof.** $X_{(k)} = x$ requires one observation at $x$ (density $f(x)$), $k-1$ observations below $x$ (each with probability $P[X \leq x] = F(x)$), and $n-k$ observations above $x$ (each with probability $P[X > x] = 1 - F(x)$). The number of such arrangements is $\dfrac{n!}{(k-1)!(n-k)!}$, since the order within the low group and within the high group is irrelevant. Multiplying gives the density. $\blacksquare$

### Corollary 3.5.1 — joint density of $X_{(i)}$ and $X_{(j)}$, $i < j$

By the same argument, for $x_{(i)} < x_{(j)}$,

$$f_{X_{(i)}, X_{(j)}}(x_{(i)}, x_{(j)}) = \frac{n!}{(i-1)!\,(j-i-1)!\,(n-j)!}[F(x_{(i)})]^{i-1}f(x_{(i)})\left[F(x_{(j)}) - F(x_{(i)})\right]^{j-i-1}f(x_{(j)})[1 - F(x_{(j)})]^{n-j}$$

### Theorem 3.5.3 — sample minimum and sample maximum

The cdfs follow directly from the definition:

$$G_1(x) = P[X_{(1)} \leq x] = 1 - P[\text{all } X_i > x] = 1 - [1 - F(x)]^n$$

$$G_n(x) = P[X_{(n)} \leq x] = P[\text{all } X_i \leq x] = [F(x)]^n$$

Differentiating in the continuous case,

$$f_{X_{(1)}}(x) = \frac{d}{dx}G_1(x) = n[1 - F(x)]^{n-1}f(x) \tag{3.13}$$

$$f_{X_{(n)}}(x) = \frac{d}{dx}G_n(x) = n[F(x)]^{n-1}f(x) \tag{3.14}$$

### Example 3.5.2 — minimum of an exponential sample

For a random sample from $f(x) = \dfrac{1}{\lambda}e^{-x/\lambda}$, $x > 0$, we have $F(x) = 1 - e^{-x/\lambda}$, so

$$f_{X_{(1)}}(x) = n\left[1 - (1 - e^{-x/\lambda})\right]^{n-1}\frac{1}{\lambda}e^{-x/\lambda} = \frac{n}{\lambda}\exp\left(-\frac{nx}{\lambda}\right), \quad x > 0 \qquad \text{(by (3.13))}$$

i.e. $X_{(1)}$ is exponential with mean $\lambda/n$.

### Definition 3.5.2 — sample median, sample range

The **sample median** is the middle order statistic, and measures location:

- $n$ odd: $\tilde{X} = X_{(k)}$, $\; k = \dfrac{n+1}{2}$
- $n$ even: $\tilde{X} = \dfrac{X_{(k)} + X_{(k+1)}}{2}$, $\; k = \dfrac{n}{2}$

The **sample range** measures dispersion:

$$R = X_{(n)} - X_{(1)}$$
