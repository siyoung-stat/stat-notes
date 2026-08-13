# 0. 기초 지식

수리통계학 본문에서 반복해서 쓰이는 미적분·대수 도구를 모은 선수 내용.

---

## 0.1 극한과 성장 속도

밀도함수의 적분이나 기댓값을 다룰 때 극한, 특히 "지수 대 다항"의 성장 속도와 로피탈 정리가 자주 쓰인다.

### 지수·로그 함수

**정의.** $e^x=\displaystyle\lim_{n\to\infty}\Big(1+\frac{x}{n}\Big)^n=\sum_{k=0}^{\infty}\frac{x^k}{k!}$, $\;(e^{ax})'=a\,e^{ax}$, $\;\displaystyle\int e^{ax}\,dx=\frac{1}{a}e^{ax}$. $\ln x$ 는 $e^x$ 의 역함수, $(\ln x)'=\tfrac1x$.

포아송 근사에서 쓰는 $\Big(1-\dfrac{\lambda}{n}\Big)^n\to e^{-\lambda}$ 도 위 극한 정의에서 나온다.

### 로피탈 정리

**정리.** $x\to a$ ($a$ 는 유한값 또는 $\pm\infty$)에서 $\dfrac{f(x)}{g(x)}$ 가 $\dfrac00$ 또는 $\dfrac{\infty}{\infty}$ 꼴 부정형이고, $a$ 근방에서 $f,g$ 가 미분가능($g'\ne0$)하며 $\displaystyle\lim\frac{f'}{g'}$ 가 존재하면

$$\lim_{x\to a}\frac{f(x)}{g(x)}=\lim_{x\to a}\frac{f'(x)}{g'(x)}.$$

**증명 (개요).** $f(a)=g(a)=0$ 이라 두자. **코시 평균값정리**에 의해 $a$ 와 $x$ 사이의 어떤 $c$ 가 존재하여

$$\frac{f(x)}{g(x)}=\frac{f(x)-f(a)}{g(x)-g(a)}=\frac{f'(c)}{g'(c)}.$$

$x\to a$ 이면 $c$ 도 $a$ 와 $x$ 사이에 끼여 $c\to a$ 이므로 $\dfrac{f(x)}{g(x)}\to\displaystyle\lim_{t\to a}\dfrac{f'(t)}{g'(t)}$. ($\infty/\infty$ 꼴도 비슷한 방식으로 증명된다.) $\blacksquare$

**예.** $\displaystyle\lim_{x\to0}\frac{\sin x}{x}\overset{0/0}{=}\lim\frac{\cos x}{1}=1$, $\;\displaystyle\lim_{x\to\infty}\frac{x}{e^x}\overset{\infty/\infty}{=}\lim\frac{1}{e^x}=0$, $\;\displaystyle\lim_{x\to0^+}x\ln x=\lim\frac{\ln x}{1/x}\overset{\infty/\infty}{=}\lim\frac{1/x}{-1/x^2}=\lim(-x)=0$.

> $0\cdot\infty,\ \infty-\infty,\ 1^\infty$ 같은 부정형은 먼저 $\dfrac00$ 또는 $\dfrac\infty\infty$ 꼴로 변형한 뒤 로피탈을 적용한다(위 $x\ln x$ 예).

### 성장 속도

**정리 (지수와 다항의 성장 비교).** 임의의 고정된 $n$ 에 대해

$$\lim_{x\to\infty}\frac{x^n}{e^x}=0,\qquad\text{즉}\quad \lim_{x\to\infty}x^n e^{-x}=0.$$

성장 속도 위계: $\ \ln x\ \ll\ x^a\ \ll\ e^x\ $ (임의의 $a\gt0$) — 지수함수가 어떤 거듭제곱·로그보다 빨리 커진다.

**증명.** $x\gt0$ 에서 $e^x=\sum_{k\ge0}\dfrac{x^k}{k!}\gt\dfrac{x^{n+1}}{(n+1)!}$ 이므로 $0\lt\dfrac{x^n}{e^x}\lt\dfrac{(n+1)!}{x}\to0$. (또는 로피탈을 $n$ 번: $\dfrac{x^n}{e^x}\to\dfrac{n!}{e^x}\to0$.) $\blacksquare$

> 이 때문에 감마함수 증명의 경계항 $\big[-t^{k-1}e^{-t}\big]_0^\infty=0$ 이 된다(지수가 다항을 누름).

---

## 0.2 적분 기법

### 미적분의 기본정리 (FTC)

$$\frac{d}{dx}\int_a^x f(t)\,dt=f(x).$$

"$a$ 부터 $x$ 까지의 넓이"를 $x$ 로 미분하면 그 끝점의 높이 $f(x)$ 가 된다. 그래서 연속형에서 cdf $F(x)=\int_{-\infty}^x f$ 를 미분하면 밀도 $f$ 가 나온다.

**증명 (개요).** $A(x)=\int_a^x f(t)\,dt$ 라 하면 $A(x+h)-A(x)=\int_x^{x+h}f(t)\,dt$. $f$ 가 연속이면 적분의 평균값정리로 이 값은 $f(\xi)\,h$ ($\xi\in[x,x+h]$). 따라서 $\dfrac{A(x+h)-A(x)}{h}=f(\xi)\to f(x)$ ($h\to0$), 곧 $A'(x)=f(x)$. $\blacksquare$

### 부분적분

**공식.**

$$\int u\,dv = uv-\int v\,du.$$

"미분하면 간단해지는 쪽"을 $u$ 로 잡는 게 요령($x^n$ → 미분하며 차수↓).

**유도.** 곱의 미분법 $(uv)'=u'v+uv'$ 의 양변을 적분하면 $uv=\int u'v\,dx+\int uv'\,dx$. 이항하면 $\int uv'\,dx=uv-\int u'v\,dx$, 즉 $\int u\,dv=uv-\int v\,du$. $\blacksquare$

**예 (감마적분과 계승).** $u=x,\ dv=e^{-x}dx$:

$$\int_0^\infty xe^{-x}dx=\big[-xe^{-x}\big]_0^\infty+\int_0^\infty e^{-x}dx=0+1=1$$

($\because$ 경계항은 성장속도로 $0$). 반복하면 $\int_0^\infty x^n e^{-x}dx=n!$. 또 $\displaystyle\int_0^\infty e^{-x/\theta}dx=\theta$, $\displaystyle\int_0^{a}ye^{-y}dy=1-(1+a)e^{-a}$ 등도 같은 방식.

### <a name="change-of-variables"></a>다변수 변수변환 (야코비안)

$(x,y)\to(u,v)$ 일대일 변환에서

$$\iint_R f(x,y)\,dx\,dy=\iint_{R'} f\big(x(u,v),y(u,v)\big)\,|J|\,du\,dv,\qquad J=\det\frac{\partial(x,y)}{\partial(u,v)}.$$

> 변환은 영역 $R$ 을 $R'$ 으로 옮긴다. 작은 넓이조각이 $|J|$ 배로 늘거나 줄기 때문에 적분에 $|J|$ 를 곱한다. 결합변환에서 밀도함수에 $|J|$ 가 곱해지는 이유다.

---

## 0.3 감마·베타·가우스 적분

감마·베타분포와 정규분포의 정규화 상수가 모두 이 절의 적분에서 나온다.

### 감마함수

**정의.**

$$\Gamma(k)=\int_0^\infty t^{k-1}e^{-t}\,dt\qquad(k\gt0).$$

**성질.**

1. $\Gamma(k)=(k-1)\,\Gamma(k-1)$
2. $\Gamma(n)=(n-1)!$ ($n$ 자연수)
3. $\Gamma(\tfrac12)=\sqrt\pi$

**증명 (1)·(2).** 부분적분 $u=t^{k-1},\ dv=e^{-t}dt$:

$$\Gamma(k)=\big[-t^{k-1}e^{-t}\big]_0^\infty+(k-1)\int_0^\infty t^{k-2}e^{-t}dt=(k-1)\Gamma(k-1)$$

(경계항은 성장속도로 $0$). $\Gamma(1)=\int_0^\infty e^{-t}dt=1$ 이므로 $\Gamma(n)=(n-1)!$. (3)은 가우스 적분에서 따라온다. $\blacksquare$

> 감마함수는 계승 $n!$ 을 실수·복소수로 확장한 것 — $\Gamma(n+1)=n!$. 감마분포 $\mathrm{GAM}(k,\theta)$ 의 정규화 상수가 $\theta^k\Gamma(k)$ 다.

### 베타함수

**정의.**

$$B(a,b)=\int_0^1 x^{a-1}(1-x)^{b-1}\,dx=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}.$$

$[0,1]$ 위 적분이라 베타분포 $\mathrm{BETA}(a,b)$ 의 정규화 상수로 쓰인다.

**유도 (베타–감마 관계).** $\Gamma(a)\Gamma(b)=\displaystyle\int_0^\infty x^{a-1}e^{-x}dx\int_0^\infty y^{b-1}e^{-y}dy=\iint_{x,y\gt0}x^{a-1}y^{b-1}e^{-(x+y)}\,dx\,dy$. 치환 $x=ts,\ y=t(1-s)$ ($t=x+y\gt0,\ s\in(0,1)$, 야코비안 $|J|=t$):

$$
\begin{aligned}
\Gamma(a)\Gamma(b)&=\int_0^\infty\!\!\int_0^1 (ts)^{a-1}\big(t(1-s)\big)^{b-1}e^{-t}\,t\,ds\,dt \\[6pt]
&=\underbrace{\int_0^\infty t^{a+b-1}e^{-t}dt}_{\Gamma(a+b)}\;\underbrace{\int_0^1 s^{a-1}(1-s)^{b-1}ds}_{B(a,b)}.
\end{aligned}
$$

따라서 $\Gamma(a)\Gamma(b)=\Gamma(a+b)\,B(a,b)$, 곧 $B(a,b)=\Gamma(a)\Gamma(b)/\Gamma(a+b)$. $\blacksquare$

### 가우스 적분

**정리 (정규분포의 정규화).**

$$\int_{-\infty}^{\infty}e^{-x^2/2}\,dx=\sqrt{2\pi}.$$

종 모양 곡선 $e^{-x^2/2}$ 아래 전체 넓이가 $\sqrt{2\pi}$. 그래서 $\tfrac{1}{\sqrt{2\pi}\sigma}e^{-(x-\mu)^2/2\sigma^2}$ 의 적분이 $1$ 이 된다.

**증명 (극좌표).**

**① 제곱해서 2차원으로.** $I=\int_{-\infty}^\infty e^{-x^2/2}dx$ 는 1차원에서 부정적분이 초등함수로 안 나온다. 대신 같은 적분을 변수 이름만 $x,y$ 로 바꿔 곱하면, 지수끼리 더해져 평면 전체의 이중적분이 된다:

$$I^2=\Big(\int_{-\infty}^\infty e^{-x^2/2}dx\Big)\Big(\int_{-\infty}^\infty e^{-y^2/2}dy\Big)=\iint_{\mathbb R^2}e^{-(x^2+y^2)/2}\,dx\,dy.$$

**② 극좌표로 바꾼다.** $x=r\cos\theta,\ y=r\sin\theta$ 이고 $x^2+y^2=r^2$. 평면 전체는 $r:0\to\infty,\ \theta:0\to2\pi$ 로 덮인다. [변수변환](#change-of-variables)에서 $dx\,dy$ 는 **야코비안 행렬식** $|J|$ 배가 되는데, $|J|$ 는 $(x,y)$ 를 새 변수 $(r,\theta)$ 로 편미분한 행렬식이다:

$$
\begin{aligned}
|J|&=\left|\det\!\begin{pmatrix}\partial x/\partial r&\partial x/\partial\theta\\[2pt]\partial y/\partial r&\partial y/\partial\theta\end{pmatrix}\right|
=\left|\det\!\begin{pmatrix}\cos\theta&-r\sin\theta\\[2pt]\sin\theta&\ \ r\cos\theta\end{pmatrix}\right| \\[6pt]
&=r\cos^2\theta+r\sin^2\theta=r.
\end{aligned}
$$

그래서 $dx\,dy=|J|\,dr\,d\theta=r\,dr\,d\theta$. (기하적으로도 반지름 $r$ 에서 각을 $d\theta$ 돌면 호 길이가 $r\,d\theta$ 라 넓이조각이 $dr\times r\,d\theta$.) 그러므로

$$I^2=\int_0^{2\pi}\!\!\int_0^\infty e^{-r^2/2}\,r\,dr\,d\theta.$$

**③ 안쪽 적분이 풀린다.** 바로 이 $r$ 덕분에 $u=r^2/2$ ($du=r\,dr$) 치환이 딱 맞아떨어진다:

$$\int_0^\infty e^{-r^2/2}\,r\,dr=\int_0^\infty e^{-u}\,du=\big[-e^{-u}\big]_0^\infty=1.$$

**④ 마무리.** 남은 $\theta$ 적분은 $\int_0^{2\pi}1\,d\theta=2\pi$. 따라서 $I^2=2\pi$, 곧 $I=\sqrt{2\pi}$. $\blacksquare$

**따름정리 (감마함수의 반정수값).** $\Gamma(\tfrac12)=\displaystyle\int_0^\infty t^{-1/2}e^{-t}dt$ 에서 $t=\tfrac{u^2}{2}$ ($dt=u\,du,\ t^{-1/2}=\sqrt2/u$)로 치환하면

$$\Gamma(\tfrac12)=\int_0^\infty\tfrac{\sqrt2}{u}e^{-u^2/2}\,u\,du=\sqrt2\int_0^\infty e^{-u^2/2}du=\sqrt2\cdot\tfrac{\sqrt{2\pi}}{2}=\sqrt\pi$$

($\because\int_0^\infty e^{-u^2/2}du=\tfrac12\sqrt{2\pi}$). $\blacksquare$

---

## 0.4 급수와 테일러 전개

### 등비급수와 지수급수

**합 공식 (등비급수).** 첫째항 $a$, 공비 $r$ 인 등비수열의 합:

1. 유한합: $\displaystyle\sum_{k=0}^{n-1}ar^k=a\,\frac{1-r^n}{1-r}\quad(r\ne1)$.
2. 무한합: $\displaystyle\sum_{k=0}^{\infty}ar^k=\frac{a}{1-r}\quad(|r|\lt1)$.
3. 미분형: $\displaystyle\sum_{k=1}^{\infty}k\,r^{k-1}=\frac{1}{(1-r)^2}\quad(|r|\lt1)$.

③은 ②($a=1$)의 양변을 $r$ 로 미분해 얻으며, 이산형 분포의 평균 계산에 자주 쓰인다.

**유도.**

**① 유한합.** $S_n=a+ar+ar^2+\cdots+ar^{n-1}$ 에 $r$ 을 곱한 $rS_n=ar+ar^2+\cdots+ar^n$ 을 빼면 가운데 항이 모두 소거되어

$$S_n-rS_n=a-ar^n=a(1-r^n)\ \Rightarrow\ S_n=a\,\frac{1-r^n}{1-r}\quad(r\ne1).$$

**② 무한합.** $|r|\lt1$ 이면 $n\to\infty$ 일 때 $r^n\to0$ 이므로

$$\sum_{k=0}^{\infty}ar^k=\lim_{n\to\infty}a\,\frac{1-r^n}{1-r}=\frac{a}{1-r}.$$

**③ 미분형.** ②에서 $a=1$ 로 둔 $\displaystyle\sum_{k=0}^{\infty}r^k=\frac{1}{1-r}$ 의 양변을 $r$ 로 미분하면 (좌변은 항별 미분)

$$\sum_{k=1}^{\infty}k\,r^{k-1}=\frac{d}{dr}\Big(\frac{1}{1-r}\Big)=\frac{1}{(1-r)^2}. \qquad \blacksquare$$

**합 공식 (지수급수).**

$$e^x=\sum_{n=0}^{\infty}\frac{x^n}{n!}\ \ (\text{모든 }x),\qquad \text{특히}\ \ \sum_{k=0}^{\infty}\frac{\lambda^k}{k!}=e^{\lambda}.$$

$e^x$ 의 테일러(매클로린) 전개다(아래 절). 포아송분포의 정규화 $\sum_x e^{-\lambda}\lambda^x/x!=1$ 이나 적률 계산에 쓰인다.

**예 (기하분포의 평균).** $X\sim\mathrm{GEO}(p)$ ($q=1-p$):

$$E(X)=\sum_{x=1}^\infty x\,p\,q^{x-1}=p\sum_{x=1}^\infty x q^{x-1}=p\cdot\dfrac{1}{(1-q)^2}=\dfrac{p}{p^2}=\dfrac1p.$$

### 테일러 급수

**정의 (테일러 급수와 테일러 다항식).** $f$ 가 $a$ 근방에서 무한번 미분가능하면, **$n$ 차 테일러 다항식** $T_n$ 과 **테일러 급수**는

$$T_n(x)=\sum_{i=0}^{n}\frac{f^{(i)}(a)}{i!}(x-a)^i,\qquad \sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(x-a)^n=f(a)+f'(a)(x-a)+\frac{f''(a)}{2!}(x-a)^2+\cdots$$

$a=0$ 이면 **매클로린 급수**다. 다항식으로 함수를 근사하는 도구이며, 차수를 올릴수록 더 넓은 범위에서 정확해진다.

**유도 (계수 $f^{(n)}(a)/n!$ 의 출처).** 테일러 급수는 "함수 $f(x)$ 가 멱급수(다항식의 무한합) 꼴로 표현된다고 가정했을 때, 계수 $c_n$ 들은 무엇이어야 하는가?"라는 질문에서 출발한다. 핵심 원리는 **점 $a$ 에서 함수 $f(x)$ 와 그 도함수들의 값을 다항식의 각 항과 일치시키는 것**이다.

**① 멱급수 표현 가정.** $f(x)$ 가 중심이 $a$ 인 멱급수로 표현된다고 가정한다:

$$f(x)=\sum_{n=0}^{\infty}c_n(x-a)^n=c_0+c_1(x-a)+c_2(x-a)^2+c_3(x-a)^3+\cdots$$

**② 양변을 차례로 미분하여 계수 $c_n$ 구하기.** 수렴반경 안에서 양변을 $x$ 에 대해 **항별 미분**한 뒤, $x=a$ 를 대입한다.

- 원래 식에 $x=a$ 대입: $\ f(a)=c_0+0+0+\cdots\ \Rightarrow\ c_0=f(a)$.
- $1$ 번 미분: $\ f'(x)=c_1+2c_2(x-a)+3c_3(x-a)^2+4c_4(x-a)^3+\cdots$ → $x=a$ 대입: $\ c_1=f'(a)$.
- $2$ 번 미분: $\ f''(x)=2\cdot1\cdot c_2+3\cdot2\cdot c_3(x-a)+4\cdot3\cdot c_4(x-a)^2+\cdots$ → $x=a$ 대입: $\ f''(a)=2!\,c_2\ \Rightarrow\ c_2=\dfrac{f''(a)}{2!}$.
- $3$ 번 미분: $\ f'''(x)=3\cdot2\cdot1\cdot c_3+4\cdot3\cdot2\cdot c_4(x-a)+\cdots$ → $x=a$ 대입: $\ f'''(a)=3!\,c_3\ \Rightarrow\ c_3=\dfrac{f'''(a)}{3!}$.

**③ 일반화.** $n$ 번 미분하면 $(x-a)^k$ 항 중 $k\lt n$ 인 항은 미분하여 사라지고, $k\gt n$ 인 항은 $x=a$ 를 넣을 때 $0$ 이 되어 다 날아간다. 오직 $k=n$ 인 항만 계수 $n!\,c_n$ 으로 남아

$$f^{(n)}(a)=n!\,c_n\ \Rightarrow\ c_n=\frac{f^{(n)}(a)}{n!}.$$

구한 계수 $c_n$ 을 ①의 멱급수에 다시 대입하면 테일러 급수 공식이 완성된다. 그래서 멱급수로 나타낼 수 있으면 계수는 이 꼴로 **유일하게** 정해진다. $\blacksquare$

> **직관 — $n$ 차 근사의 의미.**
> - $1$ 차 근사 $T_1(x)=f(a)+f'(a)(x-a)$: $a$ 에서 **함숫값**과 **기울기**를 동일하게 맞춘 직선(접선).
> - $2$ 차 근사 $T_2(x)=f(a)+f'(a)(x-a)+\dfrac{f''(a)}{2!}(x-a)^2$: 휘어짐(곡률, $2$ 차 미분)까지 똑같이 맞춘 포물선.
> - $n$ 차 근사: $n$ 차 도함수 정보까지 끌어써서 $a$ 근방에서 함수를 복원해내는 다항식.

**테일러 정리 (다항식과 나머지).** $f(x)=T_n(x)+R_n(x)$ 로 쓰고 $R_n$ 을 **나머지**(remainder)라 한다. $f$ 가 자기 테일러 급수와 같을 필요충분조건은 나머지가 사라지는 것이다:

$$f(x)=\sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(x-a)^n\ \iff\ \lim_{n\to\infty}R_n(x)=0.$$

급수가 수렴하더라도 그 합이 $f$ 라는 보장은 없으므로, **$R_n\to0$ 을 반드시 확인**해야 한다.

**성질 (나머지 항 — 라그랑주형·적분형·테일러 부등식).**

- **라그랑주형**: $a$ 와 $x$ 사이의 어떤 $z$ 에 대해 $R_n(x)=\dfrac{f^{(n+1)}(z)}{(n+1)!}(x-a)^{n+1}$ ($n=0$ 이면 평균값정리).
- **적분형**: $R_n(x)=\dfrac{1}{n!}\displaystyle\int_a^x (x-t)^n\,f^{(n+1)}(t)\,dt$.
- **테일러 부등식**: $|x-a|\le d$ 에서 $|f^{(n+1)}|\le M$ 이면 $|R_n(x)|\le\dfrac{M}{(n+1)!}\,|x-a|^{n+1}$.

**주요 매클로린 전개.**

- $e^x=\displaystyle\sum_{n\ge0}\frac{x^n}{n!}$ (모든 $x$).
- $\sin x=\displaystyle\sum_{n\ge0}\frac{(-1)^n}{(2n+1)!}\,x^{2n+1}$ (모든 $x$).
- $\cos x=\displaystyle\sum_{n\ge0}\frac{(-1)^n}{(2n)!}\,x^{2n}$ (모든 $x$).
- $\dfrac{1}{1-x}=\displaystyle\sum_{n\ge0}x^n$ ($|x|\lt1$).
- $\ln(1+x)=\displaystyle\sum_{n\ge1}\frac{(-1)^{n-1}}{n}\,x^n$ ($-1\lt x\le1$).
- $(1+x)^\alpha=\displaystyle\sum_{n\ge0}\binom{\alpha}{n}x^n$ ($|x|\lt1$).

**유도 ($e^x$·$\sin x$·$\cos x$ 전개와 수렴).**

- $e^x$: $f^{(n)}(x)=e^x$ 라 $f^{(n)}(0)=1$ → $\sum_{n\ge0}x^n/n!$. **수렴**: $|x|\le d$ 에서 $|f^{(n+1)}|\le e^d$ 이므로 테일러 부등식으로 $|R_n|\le\dfrac{e^d}{(n+1)!}|x|^{n+1}\to0$ (모든 $x$).
- $\sin x$: 도함수가 $\cos x,-\sin x,-\cos x,\sin x$ 로 4주기라 $f^{(n)}(0)=0,1,0,-1,\dots$ → 홀수 차수만 남아 $x-\dfrac{x^3}{3!}+\dfrac{x^5}{5!}-\cdots$. $|f^{(n+1)}|\le1$ 이라 $R_n\to0$ (모든 $x$).
- $\cos x$: $\sin x$ 의 급수를 **항별 미분**하면 곧바로 얻는다.

$\blacksquare$

> **적률생성함수와의 연결:** $M_X(t)=E(e^{tX})=E\Big(\sum_r\dfrac{(tX)^r}{r!}\Big)=\sum_r\dfrac{t^r}{r!}E(X^r)$. 이 테일러 전개의 $t^r$ 계수가 $\dfrac{E(X^r)}{r!}$ 라서 $M_X^{(r)}(0)=E(X^r)$ 가 나온다. 확률생성함수도 같은 원리($s^j$ 계수 $=P(X=j)$).

---

## 0.5 다변량 미분과 그라디언트

여러 변수의 함수 $f(x_1,\dots,x_n)$ 를 다루는 미분이다. 최대가능도추정(MLE)처럼 여러 모수를 동시에 최적화하거나, 결합분포·변수변환을 다룰 때 바탕이 된다.

### 편도함수

**정의.** 나머지 변수를 상수로 고정하고 한 변수 $x_i$ 로만 미분한 것:

$$\frac{\partial f}{\partial x_i}=\lim_{h\to0}\frac{f(\dots,x_i+h,\dots)-f(\dots,x_i,\dots)}{h}.$$

### 그라디언트

**정의.** 모든 편도함수를 모은 벡터 $\nabla f$:

$$\nabla f=\Big(\frac{\partial f}{\partial x_1},\ \frac{\partial f}{\partial x_2},\ \dots,\ \frac{\partial f}{\partial x_n}\Big).$$

**방향도함수**: 단위벡터 $\mathbf u$ 방향의 변화율은 $D_{\mathbf u}f=\nabla f\cdot\mathbf u$.

$\nabla f$ 는 등고선(같은 값의 선)에 **수직**이고, 함숫값이 **가장 빠르게 커지는 방향**을 가리킨다. 크기 $\|\nabla f\|$ 가 그 최대 변화율.

> **최적화 연결:** 극대·극소에서는 모든 방향으로 변화율이 $0$ 이라 $\nabla f=\mathbf 0$ (임계점). 가능도함수의 $\nabla=\mathbf0$ 을 푸는 것이 곧 MLE다.

### 헤시안과 연쇄법칙

**헤시안.** 2차 편도함수 행렬 $H=\Big[\dfrac{\partial^2 f}{\partial x_i\,\partial x_j}\Big]$. 임계점에서 $H$ 가 양의정부호면 극소, 음의정부호면 극대(2차 판정).

**다변수 연쇄법칙.** $x=x(t),\,y=y(t)$ 일 때

$$\frac{d}{dt}f(x(t),y(t))=\frac{\partial f}{\partial x}x'(t)+\frac{\partial f}{\partial y}y'(t)=\nabla f\cdot(x',y').$$

### 야코비 행렬과 부피 배율

**야코비 행렬.** 벡터함수 $\mathbf g=(g_1,\dots,g_n):\mathbb R^n\to\mathbb R^n$ 의 1차 도함수 행렬과 그 행렬식(야코비안):

$$J=\Big[\frac{\partial g_i}{\partial x_j}\Big]_{n\times n},\qquad \det J.$$

**선형대수 원리 (행렬식 = 부피 배율).** 선형사상 $A$ 는 단위정육면체를 부피 $|\det A|$ 의 평행육면체로 보낸다 — 즉 $|\det A|$ 가 **부피(넓이)가 늘어나는 배율**이다. 비선형 변환은 한 점 근방에서 야코비 행렬 $J$ 로 **국소 선형근사**되므로, 작은 부피조각이 $|\det J|$ 배로 바뀐다.

> 그래서 변수변환 적분에서 $dx\,dy=|\det J|\,du\,dv$ 가 되고([0.2절](#change-of-variables)), 밀도함수에 $|J|$ 가 곱해진다. (행렬식의 기하적 의미는 선형대수의 핵심 주제다.)

---

## 0.6 삼각부등식

"합의 크기는 크기의 합을 넘지 못한다"는 부등식이다. 절댓값·노름·적분·기댓값 어디서나 같은 꼴로 쓰이며, 기댓값의 존재($E|X|\lt\infty$)나 여러 부등식 증명의 바탕이 된다.

### 실수와 벡터

**삼각부등식.**

$$|a+b|\le|a|+|b|\qquad(\text{벡터}:\ \|\mathbf x+\mathbf y\|\le\|\mathbf x\|+\|\mathbf y\|).$$

벡터 $\mathbf x+\mathbf y$ 의 길이는 $\mathbf x$, $\mathbf y$ 를 이어 간 길이의 합 이하 — 직선이 가장 짧기 때문.

**증명 (실수).** $-|a|\le a\le|a|$, $-|b|\le b\le|b|$ 를 더하면 $-(|a|+|b|)\le a+b\le|a|+|b|$. 곧 $|a+b|\le|a|+|b|$.

(벡터는 $\|\mathbf x+\mathbf y\|^2=\|\mathbf x\|^2+2\,\mathbf x\!\cdot\!\mathbf y+\|\mathbf y\|^2\le\|\mathbf x\|^2+2\|\mathbf x\|\|\mathbf y\|+\|\mathbf y\|^2=(\|\mathbf x\|+\|\mathbf y\|)^2$, 여기서 $\mathbf x\!\cdot\!\mathbf y\le\|\mathbf x\|\|\mathbf y\|$ 는 코시–슈바르츠.) $\blacksquare$

**역삼각부등식.**

$$\big|\,|a|-|b|\,\big|\le|a-b|.$$

**증명.** $|a|=|(a-b)+b|\le|a-b|+|b|$ 이라 $|a|-|b|\le|a-b|$. $a,b$ 를 바꾸면 $|b|-|a|\le|a-b|$. 두 식을 합치면 결론. $\blacksquare$

### 적분·기댓값 형태

$$\Big|\int f(x)\,dx\Big|\le\int|f(x)|\,dx,\qquad |E(X)|\le E\big(|X|\big).$$

**증명.** 모든 $x$ 에서 $-|f(x)|\le f(x)\le|f(x)|$ 이므로 적분의 단조성으로 $-\int|f|\le\int f\le\int|f|$, 곧 $|\int f|\le\int|f|$. $f$ 를 밀도와 곱한 $x f(x)$ 에 적용하면 $|E(X)|=\big|\int x f\big|\le\int|x|f=E|X|$ — 그래서 기댓값을 정의할 때 $E|X|\lt\infty$(절대수렴)를 요구한다. $\blacksquare$

> 이 부등식은 일반화 삼각부등식 $\big|\sum_i a_i\big|\le\sum_i|a_i|$ 의 연속판이다. 부울 부등식, 마르코프 부등식 등 여러 증명에서 "절댓값을 안으로 넣는" 단계로 쓰인다.

---

## 0.7 이차식과 판별식

이차식의 근의 개수와 부호를 한 수로 판정하는 도구. "이차식이 모든 값에서 $\ge0$"을 쓰는 증명(코시–슈바르츠 등)에서 핵심이다.

**정의 (판별식).** 이차식 $at^2+bt+c$ ($a\ne0$) 의 판별식은

$$D=b^2-4ac.$$

**성질 (근의 개수).**

1. $D\gt0$: 서로 다른 두 실근
2. $D=0$: 중근(한 실근)
3. $D\lt0$: 실근 없음

$a\gt0$(아래로 볼록)일 때: $D\gt0$ 면 $x$ 축과 두 번 만나고, $D=0$ 이면 접하고, $D\lt0$ 이면 만나지 않아 항상 축 위(양수)에 있다.

**성질 (이차식의 부호).** $a\gt0$ 이면 이차식은 아래로 볼록(최솟값을 가짐)이므로

$$at^2+bt+c\ge0\ \ (\text{모든 } t)\iff D=b^2-4ac\le0.$$

($a\lt0$ 이면 위로 볼록이라 "항상 $\le0$"와 같은 조건이 된다.)

> **활용 — 코시–슈바르츠 부등식:** $0\le E[(tX-Y)^2]=t^2E(X^2)-2t\,E(XY)+E(Y^2)$ 가 **모든 $t$ 에서** $\ge0$ 이므로, $t$ 에 대한 이 이차식의 판별식이 $\le0$:
>
> $$4[E(XY)]^2-4E(X^2)E(Y^2)\le0\ \Rightarrow\ [E(XY)]^2\le E(X^2)E(Y^2).$$
