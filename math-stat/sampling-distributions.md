# 3. 표본분포와 근사

## 3.1 모집단과 표본

### <a name="def-3-1-1"></a>정의 3.1.1 — 모집단, 표본

- **모집단**: 연구 대상 전체
- **표본**: 실제로 관측된 모집단의 일부

### <a name="def-3-1-2"></a>정의 3.1.2 — 확률표본

결합밀도가 $f$ 의 곱으로 분해되면 $X_1, \dots, X_n$ 을 밀도 $f$ 인 모집단에서 뽑은 크기 $n$ 의 **확률표본**이라 한다.

$$f_{X_1,\dots,X_n}(x_1, \dots, x_n) = f(x_1)f(x_2)\cdots f(x_n) = \prod_{i=1}^{n} f(x_i)$$

즉 $X_i$ 는 **iid**(독립이고 동일한 분포)이다.

### <a name="ex-3-1-1"></a>예 3.1.1 — 베르누이 모집단

$f(x) = p^x q^{1-x}$, $x = 0, 1$, $p + q = 1$, $n = 10$:

$$f_{X_1,\dots,X_{10}}(x_1, \dots, x_{10}) = \prod_{i=1}^{10} p^{x_i} q^{1-x_i} = p^{\sum x_i} q^{10 - \sum x_i} \qquad \text{(독립성)}$$

### <a name="def-3-1-3"></a>정의 3.1.3 — 통계량

**통계량**은 미지의 모수를 포함하지 않는 확률표본의 함수 $T = T(X_1, \dots, X_n)$ 이다. 확률변수의 함수이므로 $T$ 자신도 확률변수다.

### <a name="ex-3-1-2"></a>예 3.1.2 — 통계량의 판별

미지의 모수 $\theta$ 에 의존하는 함수는 통계량이 아니다.

- 통계량: $\displaystyle \bar{X}_n = \frac{1}{n}\sum_{i=1}^{n} X_i$, $\quad \max\{X_1, \dots, X_n\}$
- 통계량 아님: $\bar{X}_n - \theta$, $\quad \max\{X_1/\theta, \dots, X_n/\theta\}$

### <a name="def-3-1-4"></a>정의 3.1.4 — 표본적률

- $r$ 차 표본적률: $\displaystyle m_r' = \frac{1}{n}\sum_{i=1}^{n} X_i^r$
- $r$ 차 중심표본적률: $\displaystyle m_r = \frac{1}{n}\sum_{i=1}^{n} (X_i - \bar{X}_n)^r$

### <a name="def-3-1-5"></a>정의 3.1.5 — 표본평균, 표본분산

- 표본평균: $\displaystyle \bar{X}_n = \frac{1}{n}\sum_{i=1}^{n} X_i \;\; (= m_1')$
- 표본분산: $\displaystyle S_n^2 = \frac{1}{n-1}\sum_{i=1}^{n} (X_i - \bar{X}_n)^2$

$m_2$ 와 달리 분모가 $n$ 이 아니라 $n-1$ 이다. 이 때문에 $E(S_n^2) = \sigma^2$ 이 되어 $S_n^2$ 이 $\sigma^2$ 의 불편추정량이 된다.

### <a name="thm-3-1-1"></a>정리 3.1.1 — $\bar{X}_n$ 과 $S_n^2$ 의 적률

평균 $\mu$, 분산 $\sigma^2$, 4차 중심적률 $\mu_4 = E(X_i - \mu)^4$ 인 모집단에서 뽑은 확률표본에 대해:

1. $E(\bar{X}_n) = \mu$
2. $\displaystyle \mathrm{Var}(\bar{X}_n) = \frac{\sigma^2}{n}$
3. $E(S_n^2) = \sigma^2$
4. $\displaystyle \mathrm{Var}(S_n^2) = \frac{1}{n}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right), \quad n > 1$

**증명.**

**(1)**

$$E(\bar{X}_n) = \frac{1}{n}\sum_{i=1}^{n} E(X_i) = \frac{1}{n} n\mu = \mu \qquad \text{(기댓값의 선형성)}$$

**(2)**

$$\mathrm{Var}(\bar{X}_n) = \frac{1}{n^2}\mathrm{Var}\left(\sum_{i=1}^{n} X_i\right) = \frac{1}{n^2}\sum_{i=1}^{n}\mathrm{Var}(X_i) = \frac{\sigma^2}{n} \qquad \text{(독립성)}$$

**(3)** 다음 항등식에서 출발한다.

$$\sum_{i=1}^{n}(X_i - \bar{X}_n)^2 = \sum_{i=1}^{n}(X_i - \mu)^2 - n(\bar{X}_n - \mu)^2$$

양변에 기댓값을 취하면

$$E\left[\sum_{i=1}^{n}(X_i - \bar{X}_n)^2\right] = n\sigma^2 - n\cdot\frac{\sigma^2}{n} = (n-1)\sigma^2$$

($\because E(X_i - \mu)^2 = \sigma^2$, $E(\bar{X}_n - \mu)^2 = \mathrm{Var}(\bar{X}_n)$)

$$\therefore\; E(S_n^2) = \frac{1}{n-1}(n-1)\sigma^2 = \sigma^2$$

**(4)** $Y_i = X_i - \mu$ 로 두면 $E(Y_i) = 0$, $E(Y_i^2) = \sigma^2$, $E(Y_i^4) = \mu_4$. $Y_i$ 가 iid이므로 서로 다른 첨자가 지수 1로 남는 항의 기댓값은 $0$ 이다.

$$W := (n-1)S_n^2 = \sum_{i=1}^{n} Y_i^2 - n\bar{Y}_n^2 \quad \text{((3)의 항등식)}, \qquad E(W) = (n-1)\sigma^2$$

$$\mathrm{Var}(S_n^2) = \frac{\mathrm{Var}(W)}{(n-1)^2} = \frac{E(W^2) - [E(W)]^2}{(n-1)^2} = \frac{E(W^2)}{(n-1)^2} - \sigma^4$$

$$W^2 = \left(\sum_{i=1}^{n} Y_i^2\right)^2 - 2n\left(\sum_{i=1}^{n} Y_i^2\right)\bar{Y}_n^2 + n^2\bar{Y}_n^4$$

**(i)** $E\left[\left(\sum_i Y_i^2\right)^2\right] = \sum_i\sum_j E(Y_i^2 Y_j^2)$ 에서 $i = j$ 인 $n$ 개 항은 $\mu_4$, $i \neq j$ 인 $n(n-1)$ 개 항은 $\sigma^4$ 를 준다(독립성):

$$E\left[\left(\sum_{i=1}^{n} Y_i^2\right)^2\right] = n\mu_4 + n(n-1)\sigma^4$$

**(ii)** $\bar{Y}_n^2 = \frac{1}{n^2}\sum_{j,k} Y_j Y_k$ 를 대입하면 $j = k$ 인 항만 남는다($\because$ 독립성과 $E(Y_j) = 0$):

$$E\left[\left(\sum_{i=1}^{n} Y_i^2\right)\bar{Y}_n^2\right] = \frac{1}{n^2}\left[n\mu_4 + n(n-1)\sigma^4\right] = \frac{\mu_4 + (n-1)\sigma^4}{n}$$

**(iii)** $\bar{Y}_n^4 = \frac{1}{n^4}\sum_{i,j,k,l} Y_i Y_j Y_k Y_l$ 에서는 두 경우만 살아남는다 — 네 첨자가 모두 같은 경우($n$ 개, $\mu_4$)와 두 쌍으로 짝지어지는 경우($3n(n-1)$ 개, $\sigma^4$):

$$E(\bar{Y}_n^4) = \frac{1}{n^4}\left[n\mu_4 + 3n(n-1)\sigma^4\right] = \frac{\mu_4 + 3(n-1)\sigma^4}{n^3}$$

(i), (ii), (iii)을 대입하면

$$E(W^2) = \left[n\mu_4 + n(n-1)\sigma^4\right] - 2n\cdot\frac{\mu_4 + (n-1)\sigma^4}{n} + n^2\cdot\frac{\mu_4 + 3(n-1)\sigma^4}{n^3}$$

$\mu_4$ 와 $\sigma^4$ 의 계수를 모으면

- $\mu_4$ 의 계수: $\displaystyle n - 2 + \frac{1}{n} = \frac{(n-1)^2}{n}$
- $\sigma^4$ 의 계수: $\displaystyle n(n-1) - 2(n-1) + \frac{3(n-1)}{n} = \frac{(n-1)(n^2 - 2n + 3)}{n}$

$$E(W^2) = \frac{(n-1)^2}{n}\mu_4 + \frac{(n-1)(n^2 - 2n + 3)}{n}\sigma^4$$

$(n-1)^2$ 로 나누고 $\sigma^4$ 를 뺀 뒤 $\sigma^4$ 항을 $n(n-1)$ 로 통분하면

$$\frac{n^2 - 2n + 3}{n(n-1)} - 1 = \frac{(n^2 - 2n + 3) - (n^2 - n)}{n(n-1)} = -\frac{n-3}{n(n-1)}$$

$$\therefore\; \mathrm{Var}(S_n^2) = \frac{E(W^2)}{(n-1)^2} - \sigma^4 = \frac{\mu_4}{n} - \frac{n-3}{n(n-1)}\sigma^4 = \frac{1}{n}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right) \qquad \blacksquare$$

---

## 3.2 정규분포와 관련된 분포

## 3.2.1 카이제곱분포

### <a name="def-3-2-1"></a>정의 3.2.1 — 카이제곱분포

모수가 $(k, \theta) = (n/2, 2)$ 인 감마분포를 자유도 $n$ 의 **카이제곱분포**라 하고 $X \sim \chi^2(n)$ 으로 쓴다.

$$f_X(x) = \frac{1}{\Gamma(n/2)\,2^{n/2}}\, x^{n/2 - 1} e^{-x/2}, \quad x > 0 \qquad (3.1)$$

### <a name="thm-3-2-1"></a>정리 3.2.1 — $\chi^2(n)$ 의 적률생성함수·평균·분산

1. $M_X(t) = (1 - 2t)^{-n/2}$, $\; t < \tfrac{1}{2}$
2. $E(X) = n$
3. $\mathrm{Var}(X) = 2n$

**증명.** $\chi^2(n) = \mathrm{GAM}(n/2, 2)$ 이고 감마분포의 적률생성함수가 $(1 - \theta t)^{-k}$ 이므로 $k = n/2$, $\theta = 2$ 를 넣으면

$$M_X(t) = (1 - 2t)^{-n/2}$$

$$E(X) = k\theta = \frac{n}{2}\cdot 2 = n, \qquad \mathrm{Var}(X) = k\theta^2 = \frac{n}{2}\cdot 2^2 = 2n \qquad \blacksquare$$

### <a name="thm-3-2-2"></a>정리 3.2.2 — 카이제곱분포의 가법성

독립인 $X_i \sim \chi^2(k_i)$, $i = 1, \dots, n$ 에 대해

$$Y = \sum_{i=1}^{n} X_i \;\sim\; \chi^2\left(\sum_{i=1}^{n} k_i\right)$$

**증명.**

$$M_Y(t) = \prod_{i=1}^{n}(1 - 2t)^{-k_i/2} = (1 - 2t)^{-\sum_i k_i/2} \qquad \text{(독립합의 mgf, 정리 3.2.1)}$$

이는 $\chi^2\left(\sum_i k_i\right)$ 의 적률생성함수이므로 적률생성함수의 유일성에 의해 결론을 얻는다. $\blacksquare$

### <a name="thm-3-2-3"></a>정리 3.2.3 — 표준정규의 제곱

$Z \sim N(0,1) \;\Rightarrow\; Y = Z^2 \sim \chi^2(1)$.

**증명.**

$$
\begin{aligned}
M_{Z^2}(t) = E(e^{tZ^2}) &= \int_{-\infty}^{\infty}\frac{1}{\sqrt{2\pi}}e^{tz^2 - z^2/2}\,dz = \int_{-\infty}^{\infty}\frac{1}{\sqrt{2\pi}}e^{-z^2(1-2t)/2}\,dz \\[6pt]
&= (1 - 2t)^{-1/2}, \quad t < \tfrac{1}{2}
\end{aligned}
$$

(가우스 적분). 이는 $\chi^2(1)$ 의 적률생성함수이므로 유일성에 의해 결론을 얻는다. $\blacksquare$

### <a name="thm-3-2-4"></a>정리 3.2.4 — 표준화된 정규의 제곱합

독립인 $X_i \sim N(\mu_i, \sigma_i^2)$, $i = 1, \dots, k$ 에 대해

$$V = \sum_{i=1}^{k}\left(\frac{X_i - \mu_i}{\sigma_i}\right)^2 \;\sim\; \chi^2(k)$$

**증명.** 표준화한 $(X_i - \mu_i)/\sigma_i \sim N(0,1)$ 은 독립이고, 각 제곱은 $\chi^2(1)$ ([정리 3.2.3](#thm-3-2-3)), 그 합은 $\chi^2(k)$ ([정리 3.2.2](#thm-3-2-2))이다. $\blacksquare$

### <a name="cor-3-2-1"></a>따름정리 3.2.1

$N(\mu, \sigma^2)$ 에서 뽑은 확률표본에 대해

$$\sum_{i=1}^{n}\frac{(X_i - \mu)^2}{\sigma^2} \;\sim\; \chi^2(n)$$

## 3.2.2 F분포

### <a name="def-3-2-2"></a>정의 3.2.2 — F분포

독립인 $U \sim \chi^2(n)$, $V \sim \chi^2(m)$ 에 대해

$$X = \frac{U/n}{V/m} \qquad (3.2)$$

의 분포를 자유도 $(n, m)$ 의 **F분포**라 하고 $X \sim F(n, m)$ 으로 쓴다. 밀도는

$$f_X(x) = \frac{\Gamma[(n+m)/2]}{\Gamma(n/2)\Gamma(m/2)}\left(\frac{n}{m}\right)^{n/2}\frac{x^{(n-2)/2}}{[1 + nx/m]^{(n+m)/2}}, \quad x > 0 \qquad (3.3)$$

**(3.3)의 유도.** 이변량 변수변환: 역변환이 $u = g(x,y)$, $v = h(x,y)$ 인 일대일 사상에 대해

$$f_{X,Y}(x,y) = f_{U,V}\big(g(x,y), h(x,y)\big)\,|J|, \qquad J = \frac{\partial(u,v)}{\partial(x,y)}$$

**(1)** $(U, V)$ 의 결합밀도는 밀도 (3.1) 두 개의 곱이다(독립성):

$$f_{U,V}(u,v) = \frac{1}{\Gamma(n/2)\Gamma(m/2)2^{(n+m)/2}}\,u^{n/2-1}v^{m/2-1}e^{-(u+v)/2}, \quad u > 0,\; v > 0$$

**(2)** 보조변수 $Y = V$ 를 붙이면 역변환과 야코비안은

$$u = \frac{n}{m}xy, \quad v = y; \qquad J = \begin{vmatrix} \dfrac{n}{m}y & \dfrac{n}{m}x \\[4pt] 0 & 1 \end{vmatrix} = \frac{n}{m}y$$

$$f_{X,Y}(x,y) = \frac{(n/m)^{n/2}}{\Gamma(n/2)\Gamma(m/2)2^{(n+m)/2}}\,x^{n/2-1}y^{(n+m)/2-1}\exp\left[-\frac{y}{2}\left(1 + \frac{nx}{m}\right)\right]$$

($y$ 의 지수: $(n/2 - 1) + (m/2 - 1) + 1 = (n+m)/2 - 1$)

**(3)** $a = \tfrac{1}{2}(1 + nx/m)$ 로 두고 $y$ 에 대해 주변화하면

$$\int_0^{\infty} y^{(n+m)/2-1}e^{-ay}\,dy = \frac{\Gamma[(n+m)/2]}{a^{(n+m)/2}} = \Gamma[(n+m)/2]\cdot\frac{2^{(n+m)/2}}{(1 + nx/m)^{(n+m)/2}}$$

(감마적분). $2^{(n+m)/2}$ 가 약분되어 (3.3)을 얻는다. $\blacksquare$

### <a name="thm-3-2-5"></a>정리 3.2.5 — $F(n, m)$ 의 평균과 분산

1. $\displaystyle E(X) = \frac{m}{m-2}$, $\; m > 2$
2. $\displaystyle \mathrm{Var}(X) = \frac{2m^2(n + m - 2)}{n(m-2)^2(m-4)}$, $\; m > 4$

**증명.** $X = \frac{m}{n}\cdot\frac{U}{V}$ 이고 $U$ 와 $1/V$ 는 독립, $m > 2$ 일 때 $E(1/V) = 1/(m-2)$ 이므로

$$E(X) = \frac{m}{n}E(U)E\!\left(\frac{1}{V}\right) = \frac{m}{n}\cdot n\cdot\frac{1}{m-2} = \frac{m}{m-2}$$

분산은 $X^2 = \dfrac{m^2}{n^2}\cdot\dfrac{U^2}{V^2}$ 에서

$$E(U^2) = \mathrm{Var}(U) + [E(U)]^2 = 2n + n^2 = n(n+2) \qquad \text{(정리 3.2.1)}$$

$$
\begin{aligned}
E\!\left(\frac{1}{V^2}\right) &= \int_0^{\infty}\frac{1}{v^2}\cdot\frac{v^{m/2-1}e^{-v/2}}{\Gamma(m/2)2^{m/2}}\,dv \\[6pt]
&= \frac{\Gamma(m/2 - 2)2^{m/2-2}}{\Gamma(m/2)2^{m/2}} = \frac{1}{(m-2)(m-4)}, \quad m > 4
\end{aligned}
$$

($\because$ 감마적분과 $\Gamma(m/2) = (m/2 - 1)(m/2 - 2)\Gamma(m/2 - 2)$)

$$E(X^2) = \frac{m^2}{n^2}E(U^2)E\!\left(\frac{1}{V^2}\right) = \frac{m^2(n+2)}{n(m-2)(m-4)} \qquad \text{(독립성)}$$

$$\mathrm{Var}(X) = E(X^2) - [E(X)]^2 = \frac{m^2(n+2)}{n(m-2)(m-4)} - \frac{m^2}{(m-2)^2} = \frac{2m^2(n + m - 2)}{n(m-2)^2(m-4)}$$

(분자: $(n+2)(m-2) - n(m-4) = 2(n + m - 2)$) $\blacksquare$

### <a name="thm-3-2-6"></a>정리 3.2.6 — F 확률변수의 역수와 백분위수 관계

$1/X = \dfrac{V/m}{U/n}$ 이므로 $X \sim F(n, m) \Rightarrow 1/X \sim F(m, n)$. 따라서 상위 백분위수에 대해

$$\frac{1}{F_{\alpha}(n, m)} = F_{1-\alpha}(m, n) \qquad (3.4)$$

이 성립하여 F분포표의 하위 부분을 따로 만들 필요가 없다.

### <a name="ex-3-2-1"></a>예 3.2.1 — F 백분위수

$X \sim F(5, 10)$ 이고 $P(X \leq a) = 0.01$ 일 때

$$a = F_{0.99}(5,10) = \frac{1}{F_{0.01}(10,5)} = \frac{1}{10.1} = 0.099 \qquad \text{((3.4)와 F분포표)}$$

## 3.2.3 t분포

### <a name="def-3-2-3"></a>정의 3.2.3 — t분포

$Z \sim N(0,1)$ 과 독립인 $U \sim \chi^2(k)$ 에 대해

$$X = \frac{Z}{\sqrt{U/k}} \qquad (3.5)$$

의 분포를 자유도 $k$ 의 **t분포**라 하고 $X \sim t(k)$ 로 쓴다. 밀도는

$$f_X(x) = \frac{\Gamma[(k+1)/2]}{\Gamma(k/2)\sqrt{k\pi}}\cdot\frac{1}{(1 + x^2/k)^{(k+1)/2}}, \quad -\infty < x < \infty \qquad (3.6)$$

**(3.6)의 유도.** 보조변수 $W = U$ 로 두고 같은 이변량 변수변환을 쓴다.

**(1)** 결합밀도(독립성):

$$f_{Z,U}(z,u) = \frac{1}{\sqrt{2\pi}}e^{-z^2/2}\cdot\frac{1}{\Gamma(k/2)2^{k/2}}u^{k/2-1}e^{-u/2}, \quad u > 0$$

**(2)** 역변환과 야코비안:

$$z = x\sqrt{\frac{w}{k}}, \quad u = w; \qquad J = \begin{vmatrix} \dfrac{\sqrt{w}}{\sqrt{k}} & \dfrac{x}{2\sqrt{kw}} \\[4pt] 0 & 1 \end{vmatrix} = \sqrt{\frac{w}{k}}$$

$$f_{X,W}(x,w) = \frac{1}{\sqrt{k\pi}\,\Gamma(k/2)2^{(k+1)/2}}\,w^{(k+1)/2-1}\exp\left[-\frac{w}{2}\left(1 + \frac{x^2}{k}\right)\right]$$

($w$ 의 지수: $w^{k/2-1}\cdot w^{1/2} = w^{(k+1)/2-1}$, 상수: $\sqrt{2\pi}\cdot 2^{k/2}\cdot\sqrt{k} = 2^{(k+1)/2}\sqrt{k\pi}$)

**(3)** $a = \tfrac{1}{2}(1 + x^2/k)$ 로 두고 $w$ 에 대해 주변화하면

$$\int_0^{\infty} w^{(k+1)/2-1}e^{-aw}\,dw = \Gamma[(k+1)/2]\cdot\frac{2^{(k+1)/2}}{(1 + x^2/k)^{(k+1)/2}}$$

$2^{(k+1)/2}$ 가 약분되어 (3.6)을 얻는다. $\blacksquare$

### <a name="ex-3-2-2"></a>예 3.2.2 — t분포의 대칭성

$X \sim t(10)$ 에서 분포표로 $P(X \leq 1.812) = 0.95$ 이고 밀도가 $x = 0$ 에 대해 대칭이므로

$$P(X \leq -1.812) = 1 - P(X \leq 1.812) = 0.05$$

### <a name="thm-3-2-7"></a>정리 3.2.7 — $t(n)$ 의 평균과 분산

1. $E(X) = 0$, $\; n > 1$
2. $\displaystyle \mathrm{Var}(X) = \frac{n}{n-2}$, $\; n > 2$

**증명.** $Z$ 가 $0$ 에 대해 대칭이고 $U$ 와 독립이므로 $X$ 도 $0$ 에 대해 대칭이다. $n > 1$ 이면 $E|X| < \infty$ 이므로 $E(X) = 0$. 그러면 $\mathrm{Var}(X) = E(X^2)$ 이고

$$X^2 = \frac{Z^2}{U/n} = nZ^2\cdot\frac{1}{U}, \qquad E(X^2) = nE(Z^2)E\!\left(\frac{1}{U}\right) \quad \text{(독립성)}$$

$$E(Z^2) = \mathrm{Var}(Z) + [E(Z)]^2 = 1$$

$$E\!\left(\frac{1}{U}\right) = \frac{1}{\Gamma(n/2)2^{n/2}}\int_0^{\infty} u^{(n/2-1)-1}e^{-u/2}\,du = \frac{\Gamma(n/2 - 1)2^{n/2-1}}{\Gamma(n/2)2^{n/2}} = \frac{1}{n-2}, \quad n > 2$$

($\because$ 감마적분 $\int_0^{\infty} u^{a-1}e^{-u/2}du = \Gamma(a)2^a$ 와 $\Gamma(n/2) = (n/2 - 1)\Gamma(n/2 - 1)$)

$$\therefore\; \mathrm{Var}(X) = E(X^2) = n\cdot 1\cdot\frac{1}{n-2} = \frac{n}{n-2} \qquad \blacksquare$$

### <a name="cor-3-2-2"></a>따름정리 3.2.2 — t 확률변수의 제곱

$X \sim t(n) \;\Rightarrow\; X^2 \sim F(1, n)$.

**증명.**

$$X^2 = \frac{Z^2}{U/n} = \frac{Z^2/1}{U/n}$$

에서 $Z^2 \sim \chi^2(1)$ ([정리 3.2.3](#thm-3-2-3))이 $U \sim \chi^2(n)$ 과 독립이므로, 이는 곧 $F(1, n)$ 의 정의 (3.2)이다. $\blacksquare$

### <a name="not-3-2-1"></a>표기 3.2.1 — 상위 백분위수

상위 꼬리확률이 $\alpha$ 인 점을 다음과 같이 쓴다.

- $z_{\alpha}$: $Z \sim N(0,1)$ 에서 $P(Z > z_{\alpha}) = \alpha$
- $t_{\alpha}(n)$: $T \sim t(n)$ 에서 $P(T > t_{\alpha}(n)) = \alpha$
- $\chi_{\alpha}^2(n)$: $X \sim \chi^2(n)$ 에서 $P(X > \chi_{\alpha}^2(n)) = \alpha$
- $F_{\alpha}(n, m)$: $F \sim F(n,m)$ 에서 $P(F > F_{\alpha}(n,m)) = \alpha$

$N(0,1)$ 과 $t(n)$ 은 $0$ 에 대해 대칭이므로

$$z_{1-\alpha} = -z_{\alpha}, \qquad t_{1-\alpha}(n) = -t_{\alpha}(n)$$

---

## 3.3 정규모집단에서의 표본추출

### <a name="thm-3-3-1"></a>정리 3.3.1 — 독립인 정규확률변수의 합

독립인 $X_i \sim N(\mu_i, \sigma_i^2)$ 에 대해

$$\sum_{i=1}^{n} X_i \;\sim\; N\left(\sum_{i=1}^{n}\mu_i,\; \sum_{i=1}^{n}\sigma_i^2\right)$$

**증명.** 정규분포의 적률생성함수가 $\exp(\mu_i t + \sigma_i^2 t^2/2)$ 이므로

$$M(t) = \prod_{i=1}^{n}\exp\left(\mu_i t + \frac{\sigma_i^2 t^2}{2}\right) = \exp\left[t\sum_{i=1}^{n}\mu_i + \frac{t^2}{2}\sum_{i=1}^{n}\sigma_i^2\right]$$

(독립합의 mgf). 이는 $N\left(\sum\mu_i, \sum\sigma_i^2\right)$ 의 적률생성함수이므로 유일성에 의해 결론을 얻는다. $\blacksquare$

### <a name="cor-3-3-1"></a>따름정리 3.3.1 — 표본평균의 분포

$N(\mu, \sigma^2)$ 에서 뽑은 크기 $n$ 의 확률표본에 대해

$$\bar{X}_n \;\sim\; N\left(\mu, \frac{\sigma^2}{n}\right)$$

([정리 3.3.1](#thm-3-3-1)과 정규분포의 선형변환 성질)

### <a name="ex-3-3-1"></a>예 3.3.1 — 생산시간

독립인 $X \sim N(6, 2^2)$, $n = 10$ 이면 $\sum X_i \sim N(60, 40)$ ([정리 3.3.1](#thm-3-3-1))이고

$$P\left(\sum_{i=1}^{10} X_i \geq 70\right) = P\left[\frac{\sum X_i - 60}{\sqrt{40}} \geq \frac{70 - 60}{\sqrt{40}}\right] = 1 - \Phi(1.58) = 0.057$$

### <a name="thm-3-3-2"></a>정리 3.3.2 — $\bar{X}_n$ 과 $S_n^2$ 의 독립성, $S_n^2$ 의 분포

$N(\mu, \sigma^2)$ 에서 뽑은 크기 $n$ 의 확률표본에 대해

1. $\bar{X}_n$ 과 $S_n^2$ 은 독립이다
2. $\displaystyle \frac{(n-1)S_n^2}{\sigma^2} \sim \chi^2(n-1)$

편차 $X_i - \bar{X}_n$ 은 제약 $\sum_i (X_i - \bar{X}_n) = 0$ 때문에 **독립이 아니므로**, $S_n^2$ 에 [정리 3.2.4](#thm-3-2-4)를 곧바로 적용할 수 없다.

**(1)의 증명 — 직교변환.** 직교행렬은 $Q^{\top}Q = QQ^{\top} = I$ 를 만족하며 제곱길이를 보존한다:

$$\|Qx\|^2 = x^{\top}Q^{\top}Qx = x^{\top}x = \|x\|^2$$

iid 표본의 공분산행렬은 $\mathrm{Cov}(X) = \sigma^2 I$ 이므로 $Y = QX$ 에 대해

$$\mathrm{Cov}(Y) = Q\,\mathrm{Cov}(X)Q^{\top} = Q(\sigma^2 I)Q^{\top} = \sigma^2 QQ^{\top} = \sigma^2 I$$

$Y$ 는 다변량정규(다변량정규의 선형변환)이고 $\Sigma = \sigma^2 I$ 이면 지수부의 이차형식이 성분별로 분리되어 결합밀도가 주변밀도의 곱이 된다. 즉 $Y_1, \dots, Y_n$ 은 서로 독립이고, 서로소인 첨자집합의 함수들도 독립이다.

첫 행이 $(1/\sqrt{n}, \dots, 1/\sqrt{n})$ 인 $Q$(**헬머트 행렬**)를 택하면

$$Y_1 = \frac{1}{\sqrt{n}}\sum_{i=1}^{n} X_i = \sqrt{n}\,\bar{X}_n$$

$$\sum_{i=2}^{n} Y_i^2 = \sum_{i=1}^{n} X_i^2 - Y_1^2 = \sum_{i=1}^{n} X_i^2 - n\bar{X}_n^2 = \sum_{i=1}^{n}(X_i - \bar{X}_n)^2 = (n-1)S_n^2$$

(길이 보존). 따라서 $\bar{X}_n = Y_1/\sqrt{n}$ 은 $Y_1$ 만의 함수이고 $S_n^2 = \frac{1}{n-1}\sum_{i=2}^{n} Y_i^2$ 은 $Y_2, \dots, Y_n$ 만의 함수인데, 첨자집합 $\{1\}$ 과 $\{2, \dots, n\}$ 이 서로소이므로 둘은 독립이다. $\blacksquare$

**(2)의 증명 — 적률생성함수.** 항등식 $\sum(X_i - \mu)^2 = \sum(X_i - \bar{X}_n)^2 + n(\bar{X}_n - \mu)^2$ 을 $\sigma^2$ 으로 나누면

$$V_1 := \sum_{i=1}^{n}\frac{(X_i - \mu)^2}{\sigma^2} = \frac{(n-1)S_n^2}{\sigma^2} + \frac{n(\bar{X}_n - \mu)^2}{\sigma^2}$$

$V_2 = \dfrac{(n-1)S_n^2}{\sigma^2}$, $V_3 = \dfrac{n(\bar{X}_n - \mu)^2}{\sigma^2}$ 로 두면 $V_1 = V_2 + V_3$ 이고

- $V_1 \sim \chi^2(n)$ ([정리 3.2.4](#thm-3-2-4))
- $V_3 = \left((\bar{X}_n - \mu)/(\sigma/\sqrt{n})\right)^2 \sim \chi^2(1)$ ([따름정리 3.3.1](#cor-3-3-1), [정리 3.2.3](#thm-3-2-3))
- $V_2$ 는 $S_n^2$ 의, $V_3$ 은 $\bar{X}_n$ 의 함수이므로 $V_2 \perp V_3$ ((1)에 의해)

$$M_{V_1}(t) = M_{V_2}(t)M_{V_3}(t) \qquad \text{(독립합의 mgf)}$$

$$M_{V_2}(t) = \frac{M_{V_1}(t)}{M_{V_3}(t)} = \frac{(1 - 2t)^{-n/2}}{(1 - 2t)^{-1/2}} = (1 - 2t)^{-(n-1)/2}, \quad t < \tfrac{1}{2} \qquad \text{(정리 3.2.1)}$$

이는 $\chi^2(n-1)$ 의 적률생성함수이므로 유일성에 의해 $V_2 \sim \chi^2(n-1)$. $\blacksquare$

회전에 대해 성분의 독립성이 유지되는 것은 정규분포뿐이므로, $\bar{X}_n$ 과 $S_n^2$ 의 독립성은 정규분포의 **특성화**다.

### <a name="ex-3-3-2"></a>예 3.3.2 — 표본분산에 대한 확률

[예 3.3.1](#ex-3-3-1)에 이어서($n = 10$, $\sigma^2 = 4$) $9S_{10}^2/4 \sim \chi^2(9)$ ([정리 3.3.2](#thm-3-3-2)(2))이므로

$$P(S_{10}^2 > 5) = P\left[\frac{9S_{10}^2}{4} > \frac{45}{4}\right] = P[\chi^2(9) > 11.25] \approx 0.259$$

### <a name="thm-3-3-3"></a>정리 3.3.3 — 두 표본분산의 비

$N(\mu_X, \sigma_X^2)$ 에서 크기 $n$, $N(\mu_Y, \sigma_Y^2)$ 에서 크기 $m$ 의 독립인 확률표본을 뽑으면

$$F = \frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2} \;\sim\; F(n-1, m-1)$$

**증명.** $(n-1)S_X^2/\sigma_X^2 \sim \chi^2(n-1)$, $(m-1)S_Y^2/\sigma_Y^2 \sim \chi^2(m-1)$ ([정리 3.3.2](#thm-3-3-2)(2))이고 두 표본이 독립이므로 이들도 독립이다. 각각을 자유도로 나누면

$$F = \frac{[(n-1)S_X^2/\sigma_X^2]/(n-1)}{[(m-1)S_Y^2/\sigma_Y^2]/(m-1)} = \frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2} \qquad \text{(정의 3.2.2)} \quad \blacksquare$$

특히 $\sigma_X^2 = \sigma_Y^2$ 이면 $S_X^2/S_Y^2 \sim F(n-1, m-1)$ 이고, 이것이 두 분산의 동일성 검정의 근거다.

### <a name="thm-3-3-4"></a>정리 3.3.4 — 스튜던트화된 t 통계량

$N(\mu, \sigma^2)$ 에서 뽑은 확률표본에 대해 $S_n = \sqrt{S_n^2}$ 이라 하면

$$T = \frac{\sqrt{n}(\bar{X}_n - \mu)}{\sqrt{\sum_{i=1}^{n}(X_i - \bar{X}_n)^2/(n-1)}} = \frac{\bar{X}_n - \mu}{S_n/\sqrt{n}} \;\sim\; t(n-1)$$

**증명.**

$$Z = \frac{\sqrt{n}(\bar{X}_n - \mu)}{\sigma} \sim N(0,1), \qquad U = \frac{(n-1)S_n^2}{\sigma^2} \sim \chi^2(n-1)$$

이고 $Z$, $U$ 는 독립이므로([정리 3.3.2](#thm-3-3-2))

$$T = \frac{Z}{\sqrt{U/(n-1)}} = \frac{\sqrt{n}(\bar{X}_n - \mu)/\sigma}{\sqrt{[(n-1)S_n^2/\sigma^2]/(n-1)}} = \frac{\sqrt{n}(\bar{X}_n - \mu)}{S_n} \;\sim\; t(n-1)$$

([정의 3.2.3](#def-3-2-3)) $\blacksquare$

**스튜던트화**는 장애모수 $\sigma^2$ 을 $S_n^2$ 으로 바꿔 제거하는 것이다. $n \to \infty$ 이면 분모가 $\sigma$ 로 수렴하므로 $t(n-1) \to N(0,1)$ (슬러츠키 정리, [예 3.4.5](#ex-3-4-5) 참고).

---

## 3.4 대수의 법칙과 중심극한정리

## 3.4.1 대수의 법칙

### <a name="def-3-4-1"></a>정의 3.4.1 — 확률수렴

$$\lim_{n\to\infty} P(|X_n - X| \geq \epsilon) = 0 \;\text{ (모든 } \epsilon > 0) \quad \Longleftrightarrow \quad X_n \xrightarrow{p} X$$

동치로 $\lim_{n\to\infty} P(|X_n - X| < \epsilon) = 1$.

### <a name="thm-3-4-1"></a>정리 3.4.1 — 대수의 법칙

평균 $\mu < \infty$ 인 모집단에서 뽑은 확률표본에 대해

$$\bar{X}_n \xrightarrow{p} \mu$$

**증명.** $\sigma^2 < \infty$ 를 가정한다. 임의의 $\epsilon > 0$ 에 대해

$$P[|\bar{X}_n - \mu| < \epsilon] \geq 1 - \frac{E(\bar{X}_n - \mu)^2}{\epsilon^2} = 1 - \frac{\sigma^2/n}{\epsilon^2} \longrightarrow 1 \quad (n \to \infty)$$

($\because$ 체비셰프 부등식과 $E(\bar{X}_n - \mu)^2 = \mathrm{Var}(\bar{X}_n) = \sigma^2/n$, [정리 3.1.1](#thm-3-1-1)). $\blacksquare$

### <a name="thm-3-4-2"></a>정리 3.4.2 — 표본분산의 확률수렴

$E(X_1^4) < \infty$ 이면

$$S_n^2 = \frac{1}{n-1}\sum_{i=1}^{n}(X_i - \bar{X}_n)^2 \xrightarrow{p} \sigma^2$$

**증명.** [정리 3.1.1](#thm-3-1-1)에서 $E(S_n^2) = \sigma^2$, $\mathrm{Var}(S_n^2) = \frac{1}{n}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right)$ 이고 $E(X_1^4) < \infty$ 는 $\mu_4 < \infty$ 를 뜻한다. 임의의 $\epsilon > 0$ 에 대해

$$P[|S_n^2 - \sigma^2| \geq \epsilon] \leq \frac{\mathrm{Var}(S_n^2)}{\epsilon^2} = \frac{1}{n\epsilon^2}\left(\mu_4 - \frac{n-3}{n-1}\sigma^4\right) \longrightarrow 0 \qquad \text{(체비셰프 부등식)}$$

괄호 안은 $\mu_4 - \sigma^4$ 로 수렴하고 $1/n$ 이 $0$ 으로 가기 때문이다. $\blacksquare$

### <a name="ex-3-4-1"></a>예 3.4.1 — 표본비율의 수렴

독립인 베르누이($p$) 확률변수에 대해 $E(X_1) = p$ 이므로

$$\hat{p}_n = \frac{1}{n}\sum_{i=1}^{n} X_i \xrightarrow{p} p \qquad \text{(정리 3.4.1)}$$

## 3.4.2 중심극한정리

### <a name="def-3-4-2"></a>정의 3.4.2 — 분포수렴

$$\lim_{n\to\infty} F_{X_n}(x) = F_X(x) \;\text{ ($F_X$ 의 모든 연속점에서)} \quad \Longleftrightarrow \quad X_n \xrightarrow{d} X$$

수렴하는 것은 확률변수가 아니라 **cdf**이며, $F_X$ 의 연속점에서만 요구한다.

### <a name="thm-3-4-3"></a>정리 3.4.3 — 연속성 정리

열린구간 $-h < t < h$ 에서 $\lim_{n\to\infty} M_n(t) = M(t)$ 이고 $M(t)$ 가 cdf $F$ 를 갖는 분포의 적률생성함수이면, $F$ 의 모든 연속점에서 $\lim_{n\to\infty} F_n(x) = F(x)$ 이다.

이산형에서는 확률생성함수를 대신 쓸 수 있다. $0 \leq s < 1$ 에서 $\lim_n G_n(s) = G(s)$ 이면 $\lim_n P(X_n = k) = a_k$ 이며, 여기서 $G(s) = \sum_k a_k s^k$.

### <a name="thm-3-4-4"></a>정리 3.4.4 — 중심극한정리

평균 $\mu$, 분산 $\sigma^2 < \infty$ 인 모집단에서 뽑은 확률표본에 대해

$$Z_n = \frac{\sum_{i=1}^{n} X_i - E\left(\sum_{i=1}^{n} X_i\right)}{\sqrt{\mathrm{Var}\left(\sum_{i=1}^{n} X_i\right)}} = \frac{\sum_{i=1}^{n}(X_i - \mu)}{\sqrt{n}\,\sigma} \xrightarrow{d} N(0,1)$$

### <a name="lem-3-4-1"></a>보조정리 3.4.1 — 1차 나머지의 라그랑주형

$f$ 가 $a$ 를 포함하는 열린구간에서 두 번 미분가능하다 하고, 나머지를 1차 근사가 놓치는 부분으로 정의하면

$$R(x) := f(x) - f(a) - f'(a)(x - a) \quad \Longleftrightarrow \quad f(x) = f(a) + f'(a)(x-a) + R(x)$$

$a$ 와 $x$ 사이의 어떤 $\xi$ 에 대해

$$R(x) = \frac{f''(\xi)}{2}(x - a)^2$$

**증명.** $x \neq a$ 를 고정하고 $R(x) = \dfrac{(x-a)^2}{2}M$ 으로 $M$ 을 정의한다.

$$g(u) = f(u) - f(a) - f'(a)(u - a) - \frac{(u-a)^2}{2}M$$

으로 두면 $g(a) = 0$, $g'(a) = 0$, $g(x) = 0$ ($M$ 의 정의에 의해)이다.

- $g(a) = g(x) = 0 \Rightarrow$ 어떤 $\xi_1 \in (a, x)$ 에 대해 $g'(\xi_1) = 0$ (롤의 정리)
- $g'(a) = g'(\xi_1) = 0 \Rightarrow$ 어떤 $\xi_2 \in (a, \xi_1)$ 에 대해 $g''(\xi_2) = 0$ (롤의 정리)

두 번 미분하면 1차 부분이 사라지고 $M$ 만 남는다:

$$g''(u) = f''(u) - M, \qquad g''(\xi_2) = f''(\xi_2) - M = 0 \iff M = f''(\xi_2)$$

$\xi = \xi_2$ 로 쓰면 결론을 얻는다. $\blacksquare$

아래에서 쓰는 요점은 $\xi$ 가 $a$ 와 $x$ 사이에 갇혀 있어 $x \to a$ 이면 $\xi \to a$ 가 강제된다는 것이다.

**[정리 3.4.4](#thm-3-4-4)의 증명 (적률생성함수의 존재를 가정).**

**(1)** 중심화한 확률변수의 적률생성함수는 모든 $i$ 에 공통이다:

$$m(t) = E\left[e^{t(X_i - \mu)}\right], \quad m(0) = 1, \quad m'(0) = E(X_i - \mu) = 0, \quad m''(0) = E(X_i - \mu)^2 = \sigma^2$$

**(2)** $0$ 근방에서 전개하고 주항을 분리하면

$$m(t) = m(0) + m'(0)t + \frac{m''(\xi)t^2}{2} = 1 + \frac{\sigma^2 t^2}{2} + \frac{[m''(\xi) - \sigma^2]t^2}{2}, \qquad 0 < \xi < t$$

(라그랑주 나머지를 쓴 테일러 전개, [보조정리 3.4.1](#lem-3-4-1))

**(3)** $t/(\sqrt{n}\sigma)$ 를 대입하면

$$
\begin{aligned}
M_{Z_n}(t) &= \prod_{i=1}^{n} M_{X_i - \mu}\left(\frac{t}{\sqrt{n}\sigma}\right) = \left[m\left(\frac{t}{\sqrt{n}\sigma}\right)\right]^n \\[6pt]
&= \left[1 + \frac{t^2}{2n} + \frac{[m''(\xi) - \sigma^2]t^2}{2n\sigma^2}\right]^n, \qquad 0 < \xi < \frac{t}{\sqrt{n}\sigma}
\end{aligned}
$$

(독립합의 mgf)

**(4)** $t$ 를 고정하면 $t/(\sqrt{n}\sigma) \to 0$ 이고 $\xi$ 가 $0$ 으로 조여진다. $m''$ 이 $0$ 에서 연속이므로 $m''(\xi) - \sigma^2 \to m''(0) - \sigma^2 = 0$. 대괄호는 $1 + (t^2/2 + o(1))/n$ 이므로

$$\lim_{n\to\infty} M_{Z_n}(t) = \exp(t^2/2) \qquad \text{(지수극한 } (1 + c/n)^n \to e^c)$$

**(5)** 이는 $N(0,1)$ 의 적률생성함수이므로 $Z_n \xrightarrow{d} N(0,1)$ ([정리 3.4.3](#thm-3-4-3)과 적률생성함수의 유일성). $\blacksquare$

중심극한정리는 모집단의 모양과 무관하게 유한한 평균과 분산만 요구한다. 평균이 없는 코시분포 등에서는 성립하지 않는다. 동치로

$$Z_n = \frac{\bar{X}_n - E(\bar{X}_n)}{\sqrt{\mathrm{Var}(\bar{X}_n)}} = \frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}} \qquad (3.8)$$

이므로 중심극한정리는 $\bar{X}_n$ 또는 $\sum X_i$ 의 분포를 정규분포로 근사한다.

### <a name="ex-3-4-2"></a>예 3.4.2 — 정규근사, 연속형

$U(0,1)$ 에서 뽑은 확률표본은 $E(X_i) = \tfrac{1}{2}$, $\mathrm{Var}(X_i) = \tfrac{1}{12}$ 이므로 $\sum X_i \approx N(n/2,\; n/12)$ 이고

$$P\left[a \leq \sum_{i=1}^{n} X_i \leq b\right] \approx \Phi\left[\frac{b - n/2}{\sqrt{n/12}}\right] - \Phi\left[\frac{a - n/2}{\sqrt{n/12}}\right] \qquad \text{(정리 3.4.4)}$$

### <a name="ex-3-4-3"></a>예 3.4.3 — 정규근사, 이항형

$X \sim B(n, p)$ 는 독립인 베르누이($p$) 확률변수의 합 $\sum_{i=1}^{n} X_i$ 와 같은 분포를 가지며 $E(X_i) = p$, $\mathrm{Var}(X_i) = pq$, $q = 1 - p$ 이다. 따라서 $\left(\sum X_i - np\right)/\sqrt{npq}$ 가 근사적으로 $N(0,1)$ 이고

$$P[a \leq X \leq b] \approx \Phi\left[\frac{b - np}{\sqrt{npq}}\right] - \Phi\left[\frac{a - np}{\sqrt{npq}}\right] \qquad \text{(정리 3.4.4)}$$

### <a name="def-3-4-3"></a>정의 3.4.3 — 연속성 수정

$P(X = k)$ 는 $k$ 를 중심으로 하는 폭 $1$ 인 막대의 넓이인데 정규근사는 곡선을 적분한다. 적분 구간을 정수 $a, b$ (막대의 중심)로 잡으면 양 끝 막대의 바깥쪽 절반이 잘린다. 각 끝을 $\tfrac{1}{2}$ 씩 늘리면 막대 전체를 덮는다:

$$P\left(a - \tfrac{1}{2} \leq X \leq b + \tfrac{1}{2}\right)$$

$$P[a \leq X \leq b] \approx \Phi\left[\frac{b + 1/2 - np}{\sqrt{npq}}\right] - \Phi\left[\frac{a - 1/2 - np}{\sqrt{npq}}\right]$$

### <a name="ex-3-4-4"></a>예 3.4.4 — 연속성 수정의 효과

$n = 50$, $p = 0.75$, $a = 36$, $b = 37$ 이면 $np = 37.5$, $npq = 9.375$.

- 수정 없이: $\Phi\left[\dfrac{37 - 37.5}{\sqrt{9.375}}\right] - \Phi\left[\dfrac{36 - 37.5}{\sqrt{9.375}}\right] = 0.1230$
- 수정 적용: $\Phi\left[\dfrac{37.5 - 37.5}{\sqrt{9.375}}\right] - \Phi\left[\dfrac{35.5 - 37.5}{\sqrt{9.375}}\right] = 0.2432$
- 이항분포 정확값: $P(36 \leq X \leq 37) = 0.2371$

## 3.4.3 슬러츠키 정리와 델타방법

### <a name="thm-3-4-5"></a>정리 3.4.5 — 슬러츠키 정리

상수 $c$ 에 대해 $X_n \xrightarrow{p} c$ 이고 $Y_n \xrightarrow{d} Z$ 이면

1. $Y_n + X_n \xrightarrow{d} Z + c$
2. $X_n Y_n \xrightarrow{d} cZ$

**(1)의 증명 (개요).** 임의의 $\epsilon > 0$ 에 대해 사건 $|X_n - c| < \epsilon$ 으로 나누면

$$F_{Y_n + X_n}(z) = P(Y_n + X_n \leq z,\; |X_n - c| < \epsilon) + P(Y_n + X_n \leq z,\; |X_n - c| \geq \epsilon)$$

$$F_{Y_n + X_n}(z) \leq P(Y_n \leq z - c + \epsilon) + P(|X_n - c| \geq \epsilon)$$

이고, 같은 분할을 $P(Y_n + c + \epsilon \leq z)$ 에 적용하면 하한을 얻으므로

$$P(Y_n \leq z - c - \epsilon) - P(|X_n - c| \geq \epsilon) \;\leq\; F_{Y_n + X_n}(z) \;\leq\; P(Y_n \leq z - c + \epsilon) + P(|X_n - c| \geq \epsilon)$$

이 부등식이 모든 $\epsilon > 0$ 에 대해 성립하고, $Z$ 가 $z - c$ 에서 연속인 곳에서는 양쪽 한계가 모두 $P(Z \leq z - c) = P(Z + c \leq z)$ 로 수렴한다($\because Y_n \xrightarrow{d} Z$, $X_n \xrightarrow{p} c$). $\blacksquare$

### <a name="ex-3-4-5"></a>예 3.4.5 — 스튜던트화 통계량의 극한분포

평균 $\mu$, 분산 $\sigma^2 < \infty$ 인 모집단에서 뽑은 확률표본(정규성은 가정하지 않음)에 대해 $S_n^2 \xrightarrow{p} \sigma^2$ ([정리 3.4.2](#thm-3-4-2))이므로 $S_n \xrightarrow{p} \sigma$ ($\because \sqrt{\cdot}$ 가 연속)이고 $(\bar{X}_n - \mu)/(\sigma/\sqrt{n}) \xrightarrow{d} N(0,1)$ ([정리 3.4.4](#thm-3-4-4))이다. 따라서

$$\frac{\bar{X}_n - \mu}{S_n/\sqrt{n}} = \frac{\sigma}{S_n}\cdot\frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}} \xrightarrow{d} N(0,1) \qquad \text{(정리 3.4.5(2))}$$

[정리 3.3.4](#thm-3-3-4)와 함께 보면 자유도가 커질수록 $t(n-1) \to N(0,1)$ 임을 확인해 준다.

### <a name="thm-3-4-6"></a>정리 3.4.6 — 델타방법

$\sqrt{n}(\bar{X}_n - \theta) \xrightarrow{d} N(0, \sigma^2)$ 이고 $g'$ 이 $\theta$ 에서 연속이며 $0$ 이 아니면

$$\sqrt{n}\left(g(\bar{X}_n) - g(\theta)\right) \xrightarrow{d} N\!\left(0,\; \sigma^2[g'(\theta)]^2\right)$$

**증명.** $\bar{X}_n$ 과 $\theta$ 사이의 어떤 $\tilde{\theta}$ 에 대해

$$
\begin{aligned}
g'(\tilde{\theta}) = \frac{g(\bar{X}_n) - g(\theta)}{\bar{X}_n - \theta} \quad &\Longleftrightarrow \quad g(\bar{X}_n) = g'(\tilde{\theta})(\bar{X}_n - \theta) + g(\theta) \\[6pt]
&\text{(평균값정리)}
\end{aligned}
$$

$$\bar{X}_n - \theta = \frac{1}{\sqrt{n}}\cdot\sqrt{n}(\bar{X}_n - \theta) \xrightarrow{d} 0\cdot N(0,\sigma^2) = 0 \qquad \text{(정리 3.4.5(2))}$$

상수로의 분포수렴은 확률수렴이므로 $\bar{X}_n \xrightarrow{p} \theta$ 이고, $\tilde{\theta}$ 가 $\bar{X}_n$ 과 $\theta$ 사이에 갇혀 있으므로 $\tilde{\theta} \xrightarrow{p} \theta$ (조임정리), $g'$ 이 연속이므로 $g'(\tilde{\theta}) \xrightarrow{p} g'(\theta)$. 따라서

$$\sqrt{n}\left(g(\bar{X}_n) - g(\theta)\right) = g'(\tilde{\theta})\cdot\sqrt{n}(\bar{X}_n - \theta) \xrightarrow{d} g'(\theta)N(0,\sigma^2) = N\!\left(0, \sigma^2[g'(\theta)]^2\right)$$

([정리 3.4.5](#thm-3-4-5)(2)) $\blacksquare$

### <a name="ex-3-4-6"></a>예 3.4.6 — 분산안정화 변환

포아송($\lambda$)에서 뽑은 확률표본은 $E(X_i) = \mathrm{Var}(X_i) = \lambda$ 이므로

$$Z_n = \sqrt{n}(\bar{X}_n - \lambda) \xrightarrow{d} N(0, \lambda) \qquad \text{(정리 3.4.4)}$$

$g(x) = \sqrt{x}$ 로 두면 $g'(\lambda) = \dfrac{1}{2\sqrt{\lambda}}$, $[g'(\lambda)]^2 = \dfrac{1}{4\lambda}$ 이므로

$$\sqrt{n}\left(\sqrt{\bar{X}_n} - \sqrt{\lambda}\right) \xrightarrow{d} N\!\left(0, \tfrac{1}{4}\right) \qquad \text{(정리 3.4.6)}$$

극한분산이 더 이상 $\lambda$ 에 의존하지 않는다. 이를 **분산안정화 변환**이라 한다.

### <a name="thm-3-4-7"></a>정리 3.4.7 — 다변량 델타방법: 근사 평균과 분산

$E(X) = \theta = (\theta_1, \dots, \theta_k)$ 인 확률벡터 $X = (X_1, \dots, X_k)$ 에 대해

$$g(x) \approx g(\theta) + \sum_{i=1}^{k} g_i'(\theta)(x_i - \theta_i), \qquad g_i'(\theta) = \left.\frac{\partial}{\partial x_i}g(x)\right|_{x = \theta} \qquad \text{(테일러 전개)}$$

$$E\,g(X) \approx g(\theta) + \sum_{i=1}^{k} g_i'(\theta)E(X_i - \theta_i) = g(\theta) \qquad (\because E(X_i - \theta_i) = 0) \qquad (3.9)$$

$$\mathrm{Var}\, g(X) \approx \sum_{i=1}^{k}[g_i'(\theta)]^2\mathrm{Var}(X_i) + 2\sum_{i > j} g_i'(\theta)g_j'(\theta)\mathrm{Cov}(X_i, X_j) \qquad (3.10)$$

(3.10)은 $\mathrm{Var}\,g(X) \approx E\left[\left(\sum_i g_i'(\theta)(X_i - \theta_i)\right)^2\right]$ 를 전개하면 나온다.

### <a name="ex-3-4-7"></a>예 3.4.7 — 비 형태의 추정량

$g(\mu_X, \mu_Y) = \dfrac{\mu_X}{\mu_Y}$ 를 $\dfrac{\bar{X}_n}{\bar{Y}_m}$ 으로 추정할 때 편도함수는

$$\frac{\partial}{\partial \mu_X}g = \frac{1}{\mu_Y}, \qquad \frac{\partial}{\partial \mu_Y}g = -\frac{\mu_X}{\mu_Y^2}$$

$$E\left(\frac{\bar{X}_n}{\bar{Y}_m}\right) \approx \frac{\mu_X}{\mu_Y} \qquad \text{((3.9))}$$

$$
\begin{aligned}
\mathrm{Var}\left(\frac{\bar{X}_n}{\bar{Y}_m}\right) &\approx \frac{1}{\mu_Y^2}\mathrm{Var}(\bar{X}_n) + \frac{\mu_X^2}{\mu_Y^4}\mathrm{Var}(\bar{Y}_m) \\[6pt]
&\quad - 2\frac{\mu_X}{\mu_Y^3}\mathrm{Cov}(\bar{X}_n, \bar{Y}_m) \qquad \text{((3.10))}
\end{aligned}
$$

$$= \left(\frac{\mu_X}{\mu_Y}\right)^2\left(\frac{\mathrm{Var}(\bar{X}_n)}{\mu_X^2} + \frac{\mathrm{Var}(\bar{Y}_m)}{\mu_Y^2} - 2\frac{\mathrm{Cov}(\bar{X}_n, \bar{Y}_m)}{\mu_X \mu_Y}\right)$$

---

## 3.5 순서통계량

### <a name="def-3-5-1"></a>정의 3.5.1 — 순서통계량

확률표본 $X_1, \dots, X_n$ 을 크기순으로 배열한 것:

$$X_{(1)} \leq X_{(2)} \leq \cdots \leq X_{(n)} \qquad (3.11)$$

### <a name="thm-3-5-1"></a>정리 3.5.1 — 순서통계량의 결합밀도

연속인 $f$ 에 대해

$$g(x_{(1)}, \dots, x_{(n)}) = \begin{cases} n!\, f(x_{(1)})f(x_{(2)})\cdots f(x_{(n)}), & x_{(1)} < x_{(2)} < \cdots < x_{(n)} \\[4pt] 0, & \text{그 외} \end{cases} \qquad (3.12)$$

**증명.** 정렬은 일대일이 아니다. $n!$ 개의 서로 다른 표본이 같은 순서벡터로 사상된다 — 예를 들어 $(5, 1, 3)$ 의 $3! = 6$ 가지 순열이 모두 $(1, 3, 5)$ 로 간다. 따라서 변수변환 공식 $f_Y(y) = f_X(g^{-1}(y))|J|$ 를 곧바로 쓸 수 없다.

표본공간을 순서가 고정된 $n!$ 개의 서로소인 영역으로 분할한다(예: $A_1 = \{x_1 < x_2 < \cdots < x_n\}$). 각 영역에서 정렬은 일대일이며 좌표를 뒤바꾸기만 하므로 $|J_i| = 1$ 이고 그 영역의 밀도는 $f(x_{(1)}) \cdots f(x_{(n)})$ 이다.

iid 표본의 밀도는 인수의 순열에 대해 불변이므로 $n!$ 개 영역이 모두 같은 기여를 한다. 이를 모두 더하면 한 영역의 밀도에 $n!$ 이 곱해져 (3.12)를 얻는다. 연속형에서는 경계 $x_i = x_j$ 의 확률이 $0$ 이므로 부등호의 강약은 밀도에 영향을 주지 않는다. $\blacksquare$

### <a name="ex-3-5-1"></a>예 3.5.1 — $f(x) = 3x^2$ 의 결합밀도

$f(x) = 3x^2$, $0 < x < 1$ 에서 뽑은 확률표본 $X_1, X_2, X_3$ 에 대해 $0 < x_{(1)} < x_{(2)} < x_{(3)} < 1$ 이면

$$g(x_{(1)}, x_{(2)}, x_{(3)}) = 3!\,(3x_{(1)}^2)(3x_{(2)}^2)(3x_{(3)}^2) = 162\,x_{(1)}^2 x_{(2)}^2 x_{(3)}^2$$

이고 그 외에서는 $0$ 이다.

### <a name="thm-3-5-2"></a>정리 3.5.2 — $k$ 번째 순서통계량의 밀도

밀도 $f$, cdf $F$ 이고 $a < x < b$ 에서 $f(x) > 0$ 인 모집단에 대해

$$f_{X_{(k)}}(x) = \frac{n!}{(k-1)!\,(n-k)!}[F(x)]^{k-1}[1 - F(x)]^{n-k}f(x), \quad a < x < b$$

이고 그 외에서는 $0$ 이다.

**증명.** $X_{(k)} = x$ 이려면 관측값 하나가 $x$ 에 있고(밀도 $f(x)$), $k-1$ 개가 $x$ 아래에 있으며(각각 확률 $P[X \leq x] = F(x)$), $n-k$ 개가 $x$ 위에 있어야 한다(각각 확률 $P[X > x] = 1 - F(x)$). 아래 집단 내부와 위 집단 내부의 순서는 무관하므로 이런 배열의 수는 $\dfrac{n!}{(k-1)!(n-k)!}$ 이다. 이들을 곱하면 밀도를 얻는다. $\blacksquare$

### <a name="cor-3-5-1"></a>따름정리 3.5.1 — $X_{(i)}$ 와 $X_{(j)}$ 의 결합밀도, $i < j$

같은 논법으로 $x_{(i)} < x_{(j)}$ 에 대해

$$
\begin{aligned}
f_{X_{(i)}, X_{(j)}}(x_{(i)}, x_{(j)}) = \frac{n!}{(i-1)!\,(j-i-1)!\,(n-j)!} &[F(x_{(i)})]^{i-1}f(x_{(i)}) \\[6pt]
&\times \left[F(x_{(j)}) - F(x_{(i)})\right]^{j-i-1}f(x_{(j)})[1 - F(x_{(j)})]^{n-j}
\end{aligned}
$$

### <a name="thm-3-5-3"></a>정리 3.5.3 — 표본최솟값과 표본최댓값

cdf는 정의에서 바로 따라온다:

$$G_1(x) = P[X_{(1)} \leq x] = 1 - P[\text{모든 } X_i > x] = 1 - [1 - F(x)]^n$$

$$G_n(x) = P[X_{(n)} \leq x] = P[\text{모든 } X_i \leq x] = [F(x)]^n$$

연속형에서 미분하면

$$f_{X_{(1)}}(x) = \frac{d}{dx}G_1(x) = n[1 - F(x)]^{n-1}f(x) \qquad (3.13)$$

$$f_{X_{(n)}}(x) = \frac{d}{dx}G_n(x) = n[F(x)]^{n-1}f(x) \qquad (3.14)$$

### <a name="ex-3-5-2"></a>예 3.5.2 — 지수분포 표본의 최솟값

$f(x) = \dfrac{1}{\lambda}e^{-x/\lambda}$, $x > 0$ 에서 뽑은 확률표본은 $F(x) = 1 - e^{-x/\lambda}$ 이므로

$$f_{X_{(1)}}(x) = n\left[1 - (1 - e^{-x/\lambda})\right]^{n-1}\frac{1}{\lambda}e^{-x/\lambda} = \frac{n}{\lambda}\exp\left(-\frac{nx}{\lambda}\right), \quad x > 0 \qquad \text{((3.13))}$$

즉 $X_{(1)}$ 은 평균 $\lambda/n$ 인 지수분포를 따른다.

### <a name="def-3-5-2"></a>정의 3.5.2 — 표본중앙값, 표본범위

**표본중앙값**은 가운데 순서통계량이며 위치를 나타낸다:

- $n$ 이 홀수: $\tilde{X} = X_{(k)}$, $\; k = \dfrac{n+1}{2}$
- $n$ 이 짝수: $\tilde{X} = \dfrac{X_{(k)} + X_{(k+1)}}{2}$, $\; k = \dfrac{n}{2}$

**표본범위**는 산포를 나타낸다:

$$R = X_{(n)} - X_{(1)}$$
