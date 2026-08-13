# 0. 기초 지식

수리통계학 본문에서 반복해 쓰이는 미적분·대수 도구.

---

## 0.1 극한과 성장 속도 <sub>limits & growth</sub>

### 지수·로그 함수 <sub>exponential & logarithm</sub>

**정의.** $e^x=\displaystyle\lim_{n\to\infty}\Big(1+\frac{x}{n}\Big)^n=\sum_{k=0}^{\infty}\frac{x^k}{k!}$, $\;(e^{ax})'=a\,e^{ax}$, $\;\displaystyle\int e^{ax}\,dx=\frac{1}{a}e^{ax}$. $\ln x$ 는 $e^x$ 의 역함수, $(\ln x)'=\tfrac1x$.

포아송 근사의 $\Big(1-\dfrac{\lambda}{n}\Big)^n\to e^{-\lambda}$ 가 위 극한 정의에서 나온다.

### 로피탈 정리 <sub>L'Hôpital's rule</sub>

**정리.** $x\to a$ ($a$ 는 유한값 또는 $\pm\infty$)에서 $\dfrac{f(x)}{g(x)}$ 가 $\dfrac00$ 또는 $\dfrac{\infty}{\infty}$ 꼴이고, $a$ 근방에서 $f,g$ 가 미분가능($g'\ne0$)하며 $\displaystyle\lim\frac{f'}{g'}$ 가 존재하면

$$\lim_{x\to a}\frac{f(x)}{g(x)}=\lim_{x\to a}\frac{f'(x)}{g'(x)}.$$

**증명 (개요).** $f(a)=g(a)=0$ 이라 두면 코시 평균값정리로 $a$ 와 $x$ 사이의 어떤 $c$ 에 대해

$$\frac{f(x)}{g(x)}=\frac{f(x)-f(a)}{g(x)-g(a)}=\frac{f'(c)}{g'(c)}.$$

$x\to a$ 이면 $c$ 가 사이에 끼여 $c\to a$. $\blacksquare$

**예.** $\displaystyle\lim_{x\to0}\frac{\sin x}{x}=1$, $\;\displaystyle\lim_{x\to\infty}\frac{x}{e^x}=0$, $\;\displaystyle\lim_{x\to0^+}x\ln x=\lim\frac{\ln x}{1/x}=\lim(-x)=0$.

$0\cdot\infty,\ \infty-\infty,\ 1^\infty$ 는 먼저 $\dfrac00$ 또는 $\dfrac\infty\infty$ 꼴로 바꾼 뒤 적용한다.

### 성장 속도 <sub>growth rates</sub>

**정리.** 임의의 고정된 $n$ 에 대해

$$\lim_{x\to\infty}\frac{x^n}{e^x}=0,\qquad \ln x\ \ll\ x^a\ \ll\ e^x\quad(a\gt0).$$

**증명.** $x\gt0$ 에서 $e^x\gt\dfrac{x^{n+1}}{(n+1)!}$ 이므로 $0\lt\dfrac{x^n}{e^x}\lt\dfrac{(n+1)!}{x}\to0$. $\blacksquare$

감마함수 증명의 경계항 $\big[-t^{k-1}e^{-t}\big]_0^\infty=0$ 이 여기서 나온다.

---

## 0.2 적분 기법 <sub>integration techniques</sub>

### 미적분의 기본정리 (FTC)

$$\frac{d}{dx}\int_a^x f(t)\,dt=f(x).$$

연속형에서 $F(x)=\int_{-\infty}^x f$ 를 미분하면 밀도 $f$ 가 나오는 근거다.

**증명 (개요).** $A(x)=\int_a^x f$ 라 하면 $A(x+h)-A(x)=\int_x^{x+h}f$. $f$ 가 연속이면 적분의 평균값정리로 이 값이 $f(\xi)h$ ($\xi\in[x,x+h]$)이므로 $\dfrac{A(x+h)-A(x)}{h}=f(\xi)\to f(x)$. $\blacksquare$

### 부분적분 <sub>integration by parts</sub>

$$\int u\,dv = uv-\int v\,du.$$

미분하면 간단해지는 쪽을 $u$ 로 잡는다.

**유도.** $(uv)'=u'v+uv'$ 의 양변을 적분해 이항. $\blacksquare$

**예.** $u=x,\ dv=e^{-x}dx$ 로 $\int_0^\infty xe^{-x}dx=\big[-xe^{-x}\big]_0^\infty+\int_0^\infty e^{-x}dx=1$ (경계항은 성장속도로 $0$). 반복하면 $\int_0^\infty x^n e^{-x}dx=n!$.

### <a name="change-of-variables"></a>다변수 변수변환 <sub>change of variables</sub>

$(x,y)\to(u,v)$ 일대일 변환에서

$$\iint_R f(x,y)\,dx\,dy=\iint_{R'} f\big(x(u,v),y(u,v)\big)\,|J|\,du\,dv,\qquad J=\det\frac{\partial(x,y)}{\partial(u,v)}.$$

넓이조각이 $|J|$ 배로 바뀌므로 적분에 $|J|$ 를 곱한다. 결합변환에서 밀도에 $|J|$ 가 곱해지는 이유.

---

## 0.3 감마·베타·가우스 적분 <sub>special integrals</sub>

감마·베타분포와 정규분포의 정규화 상수가 모두 여기서 나온다.

### 감마함수 <sub>gamma function</sub>

**정의.**

$$\Gamma(k)=\int_0^\infty t^{k-1}e^{-t}\,dt\qquad(k\gt0).$$

계승 $n!$ 을 실수로 확장한 것 — $\Gamma(n+1)=n!$.

**성질.**

1. $\Gamma(k)=(k-1)\,\Gamma(k-1)$
2. $\Gamma(n)=(n-1)!$ ($n$ 자연수)
3. $\Gamma(\tfrac12)=\sqrt\pi$

**증명 (1)·(2).** 부분적분 $u=t^{k-1},\ dv=e^{-t}dt$:

$$\Gamma(k)=\big[-t^{k-1}e^{-t}\big]_0^\infty+(k-1)\int_0^\infty t^{k-2}e^{-t}dt=(k-1)\Gamma(k-1)$$

(경계항은 $0$). $\Gamma(1)=1$ 이므로 $\Gamma(n)=(n-1)!$. (3)은 가우스 적분에서 따라온다. $\blacksquare$

### 베타함수 <sub>beta function</sub>

**정의.**

$$B(a,b)=\int_0^1 x^{a-1}(1-x)^{b-1}\,dx=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}.$$

$[0,1]$ 위 적분이라 베타분포의 정규화 상수가 된다.

**유도.** $\Gamma(a)\Gamma(b)=\iint_{x,y\gt0}x^{a-1}y^{b-1}e^{-(x+y)}\,dx\,dy$ 에 $x=ts,\ y=t(1-s)$ ($t\gt0,\ s\in(0,1)$, $|J|=t$)를 치환하면

$$\Gamma(a)\Gamma(b)=\underbrace{\int_0^\infty t^{a+b-1}e^{-t}dt}_{\Gamma(a+b)}\;\underbrace{\int_0^1 s^{a-1}(1-s)^{b-1}ds}_{B(a,b)}. \qquad \blacksquare$$

### 가우스 적분 <sub>Gaussian integral</sub>

**정리.**

$$\int_{-\infty}^{\infty}e^{-x^2/2}\,dx=\sqrt{2\pi}.$$

**증명 (극좌표).** 1차원에서는 부정적분이 초등함수로 안 나오므로 제곱해 이중적분으로 만든다:

$$I^2=\iint_{\mathbb R^2}e^{-(x^2+y^2)/2}\,dx\,dy.$$

$x=r\cos\theta,\ y=r\sin\theta$ 로 바꾸면

$$|J|=\left|\det\!\begin{pmatrix}\cos\theta&-r\sin\theta\\[2pt]\sin\theta&\ \ r\cos\theta\end{pmatrix}\right|=r,$$

곧 $dx\,dy=r\,dr\,d\theta$. 이 $r$ 덕분에 $u=r^2/2$ 치환이 맞아떨어진다:

$$I^2=\int_0^{2\pi}\!\!\int_0^\infty e^{-r^2/2}\,r\,dr\,d\theta=2\pi\int_0^\infty e^{-u}du=2\pi. \qquad \blacksquare$$

**따름정리.** $\Gamma(\tfrac12)=\int_0^\infty t^{-1/2}e^{-t}dt$ 에 $t=\tfrac{u^2}{2}$ 를 치환하면

$$\Gamma(\tfrac12)=\sqrt2\int_0^\infty e^{-u^2/2}du=\sqrt2\cdot\tfrac{\sqrt{2\pi}}{2}=\sqrt\pi. \qquad \blacksquare$$

---

## 0.4 급수와 테일러 전개 <sub>series & Taylor</sub>

### 등비급수와 지수급수 <sub>geometric & exponential series</sub>

**등비급수.** 첫째항 $a$, 공비 $r$:

1. $\displaystyle\sum_{k=0}^{n-1}ar^k=a\,\frac{1-r^n}{1-r}\quad(r\ne1)$
2. $\displaystyle\sum_{k=0}^{\infty}ar^k=\frac{a}{1-r}\quad(|r|\lt1)$
3. $\displaystyle\sum_{k=1}^{\infty}k\,r^{k-1}=\frac{1}{(1-r)^2}\quad(|r|\lt1)$

**유도.** ①은 $S_n-rS_n=a(1-r^n)$ 에서, ②는 $|r|\lt1$ 일 때 $r^n\to0$ 에서, ③은 ②($a=1$)의 양변을 $r$ 로 항별 미분해서 얻는다. $\blacksquare$

③은 이산형 분포의 평균 계산에 쓰인다. 예로 $X\sim\mathrm{GEO}(p)$ 이면

$$E(X)=p\sum_{x=1}^\infty x q^{x-1}=\frac{p}{(1-q)^2}=\frac1p.$$

**지수급수.**

$$e^x=\sum_{n=0}^{\infty}\frac{x^n}{n!}\quad(\text{모든 }x).$$

포아송분포의 정규화 $\sum_x e^{-\lambda}\lambda^x/x!=1$ 에 쓰인다.

### 테일러 급수 <sub>Taylor series</sub>

**정의.** $f$ 가 $a$ 근방에서 무한번 미분가능할 때, $n$ 차 테일러 다항식과 테일러 급수는

$$T_n(x)=\sum_{i=0}^{n}\frac{f^{(i)}(a)}{i!}(x-a)^i,\qquad \sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(x-a)^n.$$

$a=0$ 이면 매클로린 급수. $T_1$ 은 접선, $T_2$ 는 곡률까지 맞춘 포물선이다.

**유도.** $f(x)=\sum_n c_n(x-a)^n$ 을 가정하고 양변을 $n$ 번 항별 미분한 뒤 $x=a$ 를 대입하면, $k\lt n$ 인 항은 미분으로 사라지고 $k\gt n$ 인 항은 $(x-a)$ 인자로 $0$ 이 되어 $k=n$ 항만 남는다:

$$f^{(n)}(a)=n!\,c_n\ \Rightarrow\ c_n=\frac{f^{(n)}(a)}{n!}. \qquad \blacksquare$$

**테일러 정리.** $f(x)=T_n(x)+R_n(x)$ 로 쓸 때

$$f(x)=\sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(x-a)^n\ \iff\ \lim_{n\to\infty}R_n(x)=0.$$

급수가 수렴해도 그 합이 $f$ 라는 보장은 없으므로 $R_n\to0$ 을 확인해야 한다.

**나머지 항.**

- 라그랑주형: $a$ 와 $x$ 사이의 어떤 $z$ 에 대해 $R_n(x)=\dfrac{f^{(n+1)}(z)}{(n+1)!}(x-a)^{n+1}$ ($n=0$ 이면 평균값정리).
- 적분형: $R_n(x)=\dfrac{1}{n!}\displaystyle\int_a^x (x-t)^n\,f^{(n+1)}(t)\,dt$.
- 테일러 부등식: $|x-a|\le d$ 에서 $|f^{(n+1)}|\le M$ 이면 $|R_n(x)|\le\dfrac{M}{(n+1)!}|x-a|^{n+1}$.

**주요 매클로린 전개.**

- $e^x=\displaystyle\sum_{n\ge0}\frac{x^n}{n!}$, $\;\sin x=\displaystyle\sum_{n\ge0}\frac{(-1)^n}{(2n+1)!}x^{2n+1}$, $\;\cos x=\displaystyle\sum_{n\ge0}\frac{(-1)^n}{(2n)!}x^{2n}$ (모든 $x$)
- $\dfrac{1}{1-x}=\displaystyle\sum_{n\ge0}x^n$, $\;(1+x)^\alpha=\displaystyle\sum_{n\ge0}\binom{\alpha}{n}x^n$ ($|x|\lt1$)
- $\ln(1+x)=\displaystyle\sum_{n\ge1}\frac{(-1)^{n-1}}{n}x^n$ ($-1\lt x\le1$)

**수렴.** $e^x$ 는 $|x|\le d$ 에서 $|f^{(n+1)}|\le e^d$, $\sin x,\cos x$ 는 $|f^{(n+1)}|\le1$ 이므로 테일러 부등식으로 $R_n\to0$ (모든 $x$).

적률생성함수 $M_X(t)=\sum_r\dfrac{t^r}{r!}E(X^r)$ 에서 $t^r$ 계수가 $E(X^r)/r!$ 라 $M_X^{(r)}(0)=E(X^r)$ 이 나온다.

---

## 0.5 다변량 미분과 그라디언트 <sub>multivariable differentiation</sub>

### 편도함수 <sub>partial derivative</sub>

**정의.** 나머지 변수를 고정하고 한 변수로만 미분한 것:

$$\frac{\partial f}{\partial x_i}=\lim_{h\to0}\frac{f(\dots,x_i+h,\dots)-f(\dots,x_i,\dots)}{h}.$$

### 그라디언트 <sub>gradient</sub>

**정의.** 편도함수를 모은 벡터

$$\nabla f=\Big(\frac{\partial f}{\partial x_1},\ \dots,\ \frac{\partial f}{\partial x_n}\Big).$$

단위벡터 $\mathbf u$ 방향의 방향도함수(directional derivative)는 $D_{\mathbf u}f=\nabla f\cdot\mathbf u$. $\nabla f$ 는 등고선에 수직이며 함숫값이 가장 빠르게 커지는 방향을 가리킨다.

극대·극소에서는 $\nabla f=\mathbf 0$ (임계점). 가능도함수의 $\nabla=\mathbf0$ 을 푸는 것이 MLE다.

### 헤시안과 연쇄법칙 <sub>Hessian & chain rule</sub>

**헤시안.** $H=\Big[\dfrac{\partial^2 f}{\partial x_i\,\partial x_j}\Big]$. 임계점에서 양의정부호면 극소, 음의정부호면 극대.

**연쇄법칙.** $x=x(t),\,y=y(t)$ 일 때

$$\frac{d}{dt}f(x(t),y(t))=\frac{\partial f}{\partial x}x'(t)+\frac{\partial f}{\partial y}y'(t)=\nabla f\cdot(x',y').$$

### 야코비 행렬 <sub>Jacobian matrix</sub>

**정의.** $\mathbf g:\mathbb R^n\to\mathbb R^n$ 의 1차 도함수 행렬과 그 행렬식:

$$J=\Big[\frac{\partial g_i}{\partial x_j}\Big]_{n\times n},\qquad \det J.$$

$|\det A|$ 는 선형사상 $A$ 가 부피를 늘리는 배율이고, 비선형 변환은 국소적으로 $J$ 로 선형근사된다. 그래서 [변수변환](#change-of-variables)에서 $dx\,dy=|\det J|\,du\,dv$ 가 된다.

---

## 0.6 삼각부등식 <sub>triangle inequality</sub>

합의 크기는 크기의 합을 넘지 못한다. 기댓값의 존재($E|X|\lt\infty$)와 여러 부등식 증명의 바탕.

### 실수와 벡터 <sub>numbers & vectors</sub>

$$|a+b|\le|a|+|b|\qquad(\text{벡터}:\ \|\mathbf x+\mathbf y\|\le\|\mathbf x\|+\|\mathbf y\|).$$

**증명.** $-|a|\le a\le|a|$, $-|b|\le b\le|b|$ 를 더하면 $-(|a|+|b|)\le a+b\le|a|+|b|$. (벡터는 $\|\mathbf x+\mathbf y\|^2\le(\|\mathbf x\|+\|\mathbf y\|)^2$ 를 코시–슈바르츠 $\mathbf x\cdot\mathbf y\le\|\mathbf x\|\|\mathbf y\|$ 로 보인다.) $\blacksquare$

**역삼각부등식**(reverse triangle inequality)**.**

$$\big|\,|a|-|b|\,\big|\le|a-b|.$$

**증명.** $|a|=|(a-b)+b|\le|a-b|+|b|$ 에서 $|a|-|b|\le|a-b|$. $a,b$ 를 바꿔 합치면 결론. $\blacksquare$

### 적분·기댓값 형태 <sub>integral & expectation forms</sub>

$$\Big|\int f(x)\,dx\Big|\le\int|f(x)|\,dx,\qquad |E(X)|\le E\big(|X|\big).$$

**증명.** $-|f|\le f\le|f|$ 에 적분의 단조성을 적용. $xf(x)$ 에 쓰면 $|E(X)|\le E|X|$ — 기댓값 정의에서 $E|X|\lt\infty$ 를 요구하는 이유다. $\blacksquare$

---

## 0.7 이차식과 판별식 <sub>quadratic & discriminant</sub>

**정의.** 이차식 $at^2+bt+c$ ($a\ne0$) 의 판별식은 $D=b^2-4ac$.

**성질.** $D\gt0$ 이면 서로 다른 두 실근, $D=0$ 이면 중근, $D\lt0$ 이면 실근이 없다. 특히 $a\gt0$ 이면

$$at^2+bt+c\ge0\ \ (\text{모든 } t)\iff D\le0.$$

**활용 — 코시–슈바르츠**(Cauchy–Schwarz)**.** $0\le E[(tX-Y)^2]=t^2E(X^2)-2t\,E(XY)+E(Y^2)$ 가 모든 $t$ 에서 성립하므로 판별식이 $\le0$:

$$[E(XY)]^2\le E(X^2)E(Y^2).$$
