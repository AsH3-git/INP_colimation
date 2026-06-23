#set text(
  font: "New Computer Modern",
  size: 14pt
)

#set par(leading: 0.5em)

#set page(
  paper: "a4",
  margin: (x: 2cm, y: 2cm),
  numbering: "1",
  number-align: center
)

#let Si = math.op("Si")
#let E = math.op("E")
#let E_1 = math.op("E_1")
#let Res = math.op("Res")

#v(10em)
#align(center)[= Вычисление дифференциального углового сечения рассеяния фотонов ОКР]
#v(2em)
#align(center)[== $"AsH"_3$ в рамках разработок для ИКИ НЦФМ]
#pagebreak()

$
omega_"max" = (E_0 lambda)/(1 + lambda), #h(1em) lambda = (4 E_0 omega_0 sin^2 alpha/2)/m^2_e, #h(1em) alpha = pi 
$
$
f(lambda, y) = (2 pi r^2_e)/(E_0 lambda)[1/(1 - y) + 1 - y - (4 y)/(lambda (1 - y)) + (4 y^2)/(lambda^2 (1 - y)^2)] = \
= (d sigma_c)/(d omega) - "сечение Клейна-Нишины"
$
$
(d sigma_c)/(d Omega) = omega_"max"/(E_0 pi) f(lambda, y(theta_gamma))/(theta^2_c (1 + (theta_gamma/theta_c)^2)^2)
$
$
"соответственно" y(theta_gamma) = omega_"max"/(E_0(1 + (theta_gamma/theta_c)^2))
$

Возникает задача вычислить полное угловое сечение фотонов ОКР для дальнейших расчетов и моделирования.

$
d sigma_c = (2 omega_"max" r^2_e)/(E^2_0 lambda theta^2_c (1 + (theta_gamma/theta_c)^2)^2) 
[1/(1 - y(theta_gamma)) + 1 - y(theta_gamma) - (4 y(theta_gamma))/(lambda (1 - y(theta_gamma))) + (4 y^2(theta_gamma))
/(lambda^2 (1 - y(theta_gamma))^2)] \
"пусть" zeta = theta^2_c (1 + (theta_gamma/theta_c)^2)^2 , "тогда, выражая через" y(theta_gamma) \
zeta = ((omega_"max" theta_c)/(y(theta_gamma) E_0))^2
$
Интеграл тогда представим полностью через переменные $y(theta_gamma)$, которую для простоты записи положим равную $y$.
$
sigma_c = limits(integral)_Omega (2 r^2_e y^2)/(lambda omega_"max" theta^2_c)[1/(1 - y) + 1 - y - (4 y)/(lambda (1 - y))
+ (4 y^2)/(lambda^2 (1 - y)^2)] d Omega = \ {(2 r^2_e)/(lambda omega_"max" theta^2_c) = A, #h(1em) 
limits(integral)^(2 pi)_0 d phi = 1} \
= A limits(integral)^(10/gamma = Gamma)_0 y^2 sin theta_gamma [1/(1 - y) + 1 - y - (4 y)/(lambda (1 - y))
+ (4 y^2)/(lambda^2 (1 - y)^2)] d theta_gamma = 
$
$
{theta_gamma = theta, #h(0.5em) theta = theta_c sqrt(omega_"max"/(E_0 y) - 1) arrow.r.double d theta = 
(- theta_c omega_"max")/(2 sqrt(omega_"max"/(E_0 y) - 1) E_0 y^2) d y} \
= - A limits(integral)^Gamma_0 (sin(theta_c sqrt(omega_"max"/(E_0 y) - 1)) theta_c omega_"max")/
(2 E_0 sqrt(omega_"max"/(E_0 y) - 1)) [1/(1 - y) + 1 - y - (4 y)/(lambda (1 - y))
+ (4 y^2)/(lambda^2 (1 - y)^2)] d y = \
{B = A (theta_c omega_"max")/E_0 = (2 r^2_e)/(lambda E_0 theta_c)} \
#v(2em)"возникли также трудности с вычислением верхнего предела" \ "по переменной игрек, поэтому обозначим его за" Delta \
= - B limits(integral)^Delta_(-infinity) sin(theta_c sqrt(omega_"max"/(E_0 y) - 1))/(2 sqrt(omega_"max"/(E_0 y) - 1))
[1/(1 - y) + 1 - y - (4 y)/(lambda (1 - y)) + (4 y^2)/(lambda^2 (1 - y)^2)] d y = \
= - B {I_1 + I_2 - I_3 - I_4 + I_5}
$
Теперь отдельно рассмотрим отдельно интегралы $I_i, i in {1, 2, 3, 4, 5}$. 
$
I_1 = limits(integral)^Delta_(-infinity) sin(theta_c sqrt(omega_"max"/(E_0 y) - 1))/
(2 (1 - y) sqrt(omega_"max"/(E_0 y) - 1)) d y = \
{omega_"max"/(E_0 y) - 1 = xi arrow.r.double y = omega_"max"/(E_0 (xi + 1)) arrow.r.double d y =
-(E_0 omega_"max")/(E^2_0 (xi + 1)^2) d xi} \
= - limits(integral)^Delta_(-infinity) (sin(theta_c sqrt(xi)) E^2_0 (xi + 1) omega_"max")/(2 (1 - omega_"max"/(E_0 (xi + 1)))
E^2_0 (xi + 1)^2 sqrt(xi)) d xi = - limits(integral)^Delta_(-infinity) (sin(theta_c sqrt(xi)) omega_"max" E_0)/(2 sqrt(xi)
(E_0 (xi + 1) - omega_"max")) d xi = \
{sqrt(xi) = psi arrow.r.double xi = psi^2 arrow.r.double d xi = 2 psi d psi} \
= - limits(integral)^Delta_(-infinity) (sin(theta_c psi) omega_"max")/(2 E_0 psi (psi^2 + 1) - 2 omega_"max" psi) d psi
$
Данный интеграл является несобственным, имеет особенность на нижнем переделе. Можно доказать, что данный интеграл сходится,
я это сделал, но сейчас писать не буду, потому что уже 3 часа ночи.

Приняв во внимание сходимость данного интеграла далее стоит задача его вычисления.
$
I_1 = -(omega_"max")/(2 E_0) limits(integral)^Delta_(-infinity) sin(theta_c psi)/(psi(psi^2 + 1 - omega_"max")) d psi 
$
Рассмотрим подынтегральную функцию в виде 
$
sin(theta_c psi) 1/a 1/(psi (psi^2 + c^2))
$
где $a = 1, b = omega_"max", c = sqrt(1 - b/a)$.
Тогда в данных обозначениях легко разложить дробь на простейшие:
$
1/(psi (psi^2 + c^2)) = 1/(psi c^2) - psi/(c^2 (psi^2 + c^2)) arrow.r.double 
I_1 = -omega_"max"/(2 E_0) {limits(integral)^Delta_(-infinity) sin(theta_c psi)/psi d psi + 
limits(integral)^Delta_(-infinity) (sin(theta_c psi) psi)/(psi^2 + c^2) d psi} \
$
Рассмотрим первый и второй интегралы по-отдельности.
$
J_1 = limits(integral)^Delta_(-infinity) sin(theta_c psi)/psi d psi = limits(integral)^0_(-infinity)
+ limits(integral)^Delta_0 = {-psi = u} = limits(integral)^(infinity)_0 sin(k u)/u d u + limits(integral)^Delta_0 
sin( k psi)/psi d psi = \ limits(=) pi/2 + Si(theta_c Delta)
$
Второй интеграл берется не так компактно
$
J_2 = limits(integral)^Delta_(-infinity) (Im(e^(i theta_c psi)) psi)/(psi^2 + c^2) d psi = 
Im limits(integral)^Delta_(-infinity) (e^(i theta_c psi) psi)/(psi^2 + c^2) d psi = K \
K = limits(integral)^infinity_(-infinity) - limits(integral)^infinity_Delta = "In"_1 -"In"_2 \
"In"_1 = limits(integral)^infinity_(-infinity) (e^(i theta_c psi) psi)/ (psi^2 + c^2) d psi
$

Контур - верхняя полуокружность. $z_0 = i c$ -- полюс первого порядка
$
limits(Res(f(psi)))_(psi = z_0) = phi(z_0)/(zeta'(z_0)) = (e^(i theta_c psi) psi)/(2)lr(|, size:
#300%)_(z=z_0) = e^(- theta_c c)/2 \
"In"_1 = 2 pi i limits(sum)_i (limits(Res(f(psi)))_(z = z_i)) = pi i e^(- theta_c c) \
"In"_2 = "по разложению дроби на простейшие над" CC = \ = 1/2(limits(integral)^infinity_Delta e^(i theta_c psi)/(psi + i c) d psi
+ limits(integral)^infinity_Delta e^(i theta_c psi)/(psi -i c) d psi)
$
Рассмотрим общий вид интегралов суммы
$
limits(integral)^infinity_Delta e^(i theta_c psi)/(psi - a) d psi = {t = -i theta_c (psi -a) arrow.r.double d psi
= i/theta_c d t} = \ = limits(integral)^infinity_Delta e^(i theta_c a)/(t e^t) d t =e^(i theta_c a) 
limits(integral)^infinity_Delta e^(-t)/t d t = e^(i theta_c a) E_1(Delta)
$
Здесь $E_1$ -- модифицированная интегральная экспонента.
$
"In"_2 = e^(- i theta_c c) E_1(z_1) + e^(i theta_c c) E_1(z_2), \
cases(
  z_1 = theta_c c - i theta_c Delta,
  z_2 -theta_c c - i theta_c Delta
)
$
$
K = i pi e^(- theta_c c) - 1/2 {e^(-i theta_c c) E_1(theta_c c - i theta_c Delta) + 
e^(i theta_c c) E_1(- theta_c c - i theta_c Delta)} = \
= (pi e^(-theta_c c))/2 - Im[e^(i theta_c Delta) E_1(i theta_c (Delta + i c))] arrow.r.double \
arrow.r.double I_1 = - (omega_"max")/(2 E_0) {pi/2 (1 + e^((-theta_c omega_"max")/(2 E_0))) Si(theta_c Delta) 
- Im[e^(i theta_c Delta) E_1(i theta_c (Delta + (i omega_"max")/(2 E_0)))]}
$
Далее предстоит взять интегралы $I_i, i in {2,3,4,5}$. Сделать это следует менее подробно, однако
не упуская ключевых моментов.
$
I_2 = limits(integral)^Delta_(-infinity) sin(theta_c sqrt(omega_"max"/(E_0 y) - 1))/(2 sqrt(omega_"max"/(E_0 y) - 1))
d y 
$
Здесь и далее первым делом будет производиться одинаковая, классическая для данного документа, замена:
$
cases(
  sqrt(omega_"max"/(E_0 y) - 1) = zeta,
  y = omega_"max"/(E_0 (zeta^2 + 1)),
  d y = -(2 omega_"max" zeta)/(E_0 (zeta^2 + 1)^2) d zeta
)
$
Интеграл тогда принимает вид:
$
- omega_"max"/E_0 limits(integral)^Delta_(-infinity) sin(theta_c zeta)/(zeta^2 + 1)^2 d zeta 
$
По разложению на простейшие дроби над $CC$ можно получить итоговый для непосредственного интегрирования вид:
$
I_2 = -omega_"max"/E_0 {i/4 limits(integral)^Delta_(-infinity) sin(theta_c zeta)/(zeta - i) d zeta - 1/4 
limits(integral)^Delta_(-infinity) sin(theta_c zeta)/(zeta - i)^2 d zeta - i/4 limits(integral)^Delta_(-infinity) 
sin(theta_c zeta)/(zeta + i) d zeta - \
- 1/4 limits(integral)^Delta_(-infinity) sin(theta_c zeta)/(zeta + i)^2 d zeta}
$
Далее введем замену $phi = zeta - i, psi = zeta + i$.
$
I_2 = -omega_"max"/E_0 {(i e^(- theta_c))/4 Im limits(integral)^Delta_(-infinity) e^(i theta_c phi)/phi d phi
- e^(- theta_c)/4 Im limits(integral)^Delta_(-infinity) e^(i theta_c phi)/phi^2 d phi - 
(i e^(theta_c))/4 Im limits(integral)^Delta_(-infinity) e^(i theta_c psi)/psi d psi - \
- e^(theta_c)/4 Im limits(integral)^Delta_(-infinity) e^(i theta_c psi)/psi^2 d psi} =
-omega_"max"/E_0 {(i e^(-theta_c))/4 Im[E(i theta_c Delta)] - (i e^(theta_c))/4 Im[E(i theta_c Delta)]
- e^(theta_c)/4 Im limits(integral)^Delta_(-infinity) e^(i theta_c psi)/psi^2 d psi - \
- e^(-theta_c)/4 Im limits(integral)^Delta_(-infinity) e^(i theta_c phi)/phi^2 d phi}
$
Интегралы со второй степенью переменных в знаменателе легко взять по частям. В общем виде это выглядит так:
$
limits(integral)^Delta_(-infinity) e^(i theta_c x)/x^2 d x = -e^(i theta_c x)/x lr(|, size:
#300%)_(-infinity)^Delta + i theta_c limits(integral)^Delta_(-infinity) e^(i theta_c x)/x d x =
-e^(i theta_c Delta)/Delta + i theta_c E(i theta_c Delta)
$
Тогда $I_2$ окончательно выглядит следующим образом:
$
I_2 = -omega_"max"/E_0 {(i e^(-theta_c))/4 Im[E(i theta_c Delta)] - (i e^(theta_c))/4 Im[E(i theta_c Delta)] 
+ e^(i theta_c Delta - theta_c)/(4 Delta) - \ 
- e^(-theta_c)/4 Re[theta_c E(i theta_c Delta)] + e^(i theta_c Delta + theta_c)/(4 Delta) - e^(theta_c)/4 
Re[theta_c E(i theta_c Delta)]}
$
$
I_3 = limits(integral)^Delta_(-infinity) (y sin(theta_c sqrt(omega_"max"/(E_0 y) - 1)))/(2 sqrt(omega_"max"/(E_0 y) - 1)) d y
= omega_"max"^2/E^2_0 limits(integral)^Delta_(-infinity) sin(theta_c zeta)/(zeta^2 + 1)^3 d zeta
$
Что по аналогии с предыдущим интегралом легко превращается в итоговое выражение для $I_3$:
$
I_3 = (omega_"max"/E_0)^2 {(Delta (3 Delta^2 + 5) sin(theta_c Delta))/(8 (Delta^2 + 1)^2) + 
(theta_c cos(theta_c Delta))/(8 (Delta^2 + 1)) - ((theta^2_c - 3 theta_c + 3) e^(theta_c))/(16) dot \
dot Re[E_1(theta_c - i theta_c Delta)] + ((theta_c^2 + 3 theta_c + 3) e^(-theta_c))/(16) Re[E_1(-theta_c - i 
theta_c Delta)]}
$
$
I_4 = limits(integral)^Delta_(-infinity) (4 y sin(theta_c sqrt(omega_"max"/(E_0 y) - 1)))/(2 lambda (1 - y) 
sqrt(omega_"max"/(E_0 y) - 1)) d y = -(4 omega^2_"max")/(E^2_0 lambda) limits(integral)^Delta_(-infinity)
sin(theta_c zeta)/(E_0 (zeta^2 + 1)^4 - omega_"max" (zeta^2 + 1)^3) d zeta
$
Раскладывая на простейшие дроби получается выражение, где все интегралы уже вычислены ранее. Ниже 
будет приведено разложение на простейшие дроби, после чего сразу будет написан ответ.
$
1/(E_0 (zeta^2 + 1)^4 - omega_"max" (zeta^2 + 1)^3) = -E^2_0/omega^3_"max" 1/(zeta^2 + 1) - 
E_0/omega^2_"max" 1/(zeta^2 + 1)^2 - 1/omega_"max" 1/(zeta^2 + 1)^3 + \
+ (E_0/omega_"max")^3 1/(zeta^2 + c^2)
$
$
I_4 = -(4 omega^2_"max")/(E^2_0 lambda) lr(\{, size: #300%) -E^2_0/omega^3_"max" (pi/2 e^(-theta_c) - 1/2[e^(-theta_c) 
Re[E_1(theta_c - i theta_c Delta)] + e^(theta_c) Re[E_1(-theta_c - i theta_c Delta)]]) - \
- E_0/omega^2_"max" ((Delta sin(theta_c Delta))/(2 (Delta^2 + 1)) + (theta_c - 1)/4 e^(theta_c) 
Re[E_1(theta_c - i theta_c Delta)] + (theta_c + 1)/4 e^(-theta_c) dot \
dot Re[E_1(-theta_c - i theta_c Delta)]) - 1/omega_"max" ((Delta (3 Delta^2 + 5) sin(theta_c Delta))/(8 
(Delta^2 + 1)^2) + (theta_c cos(theta_c Delta))/(8 (Delta^2 + 1)) - ((theta^2_c - 3 theta_c + 3)e^(theta_c))/16 dot
$
$
dot Re[E_1(theta_c - i theta_c Delta)] + ((theta^2_c + 3 theta_c + 3) e^(-theta_c))/(16) 
Re[E_1(- theta_c - i theta_c Delta)]) + E^3_0/omega^3_"max" dot \ 
dot (pi/(2 c) e^(-theta_c c) - 1/(2 c) [e^(-theta_c c) Re[E_1(theta_c c - i theta_c Delta)] + e^(theta_c c) 
Re[E_1(-theta_c c - i theta_c Delta)]])lr(}, size: #300%)
$
$
I_5 = limits(integral)^Delta_(-infinity) sin(theta_c sqrt(omega_"max"/(E_0 y) - 1))/(2 
sqrt(omega_"max"/(E_0 y) - 1)) (4 y^2)/(lambda^2 (1 - y)^2) d y = 
-(4 omega^3_"max")/(E^2_0 lambda^2) limits(integral)^Delta_(-infinity) sin(theta_c zeta)/(E_0 (zeta^2 + 1)^4 
-omega_"max" (zeta^2 + 1)^2) d zeta
$
Разложение подынтегральной функции без учета синуса в числителе выглядит следующим образом:
$
1/D(zeta) = -1/omega_"max" 1/(zeta^2 + 1)^2 + sqrt(E_0)/(2 omega^(3/2)_"max") 1/(zeta^2 + a^2)
- sqrt(E_0)/(2 omega^(3/2)_"max") 1/(zeta^2 + b^2)
$
Здесь используются новые параметры
$
cases(
  r = sqrt(omega_"max"/E_0),
  a = sqrt(1 - r),
  b = sqrt(1 + r)
)
$
Тогда пятый интеграл записывается так:
$
I_5 = -(4 omega^3_"max")/(E^2_0 lambda^2) {-1/omega_"max" ((Delta sin(theta_c Delta))/(2(Delta^2 + 1)) + 
(theta_c cos(theta_c Delta))/(8(Delta^2 + 1)) - ((theta_c^2 - 3theta_c + 3)e^(theta_c))/(16) dot \
dot Re[E_1(theta_c - i theta_c Delta)] + ((theta^2_c - 3theta_c + 3)e^(-theta_c))/(16) 
Re[E_1(-theta_c - i theta_c Delta)]) + sqrt(E_0)/(omega^(3/2)_"max") dot \
dot (pi/(2 a) e^(-theta_c a) -1/(2a)[e^(-theta_c a) Re[E_1(theta_c a - i theta_c Delta)]
+ e^(theta_c a) Re[E_1(-theta_c a) - i theta_c Delta]] + \
+ pi/(2 b) e^(-theta_c b) - 1/(2 b) [e^(-theta_c b) Re[E_1(theta_c b - i theta_c Delta)] + e^(theta_c b) 
Re[E_1(- theta_c b - i theta_c Delta)]])}
$
Тогда итоговый интеграл, определяющий полное угловое сечение фотонов ОКР определяется следующим выражением:
$
sigma_Omega = -(2 r^2_e)/(lambda E_0 theta_c){- (omega_"max")/(2 E_0) {pi/2 (1 + e^((-theta_c omega_"max")/(2 E_0))) Si(theta_c Delta) 
- Im[e^(i theta_c Delta) E_1(i theta_c (Delta + (i omega_"max")/(2 E_0)))]} -\
-omega_"max"/E_0 {(i e^(-theta_c))/4 Im[E(i theta_c Delta)] - (i e^(theta_c))/4 Im[E(i theta_c Delta)] 
+ e^(i theta_c Delta - theta_c)/(4 Delta) - \ 
- e^(-theta_c)/4 Re[theta_c E(i theta_c Delta)] + e^(i theta_c Delta + theta_c)/(4 Delta) - e^(theta_c)/4 
Re[theta_c E(i theta_c Delta)]} -\
- (omega_"max"/E_0)^2 {(Delta (3 Delta^2 + 5) sin(theta_c Delta))/(8 (Delta^2 + 1)^2) + 
(theta_c cos(theta_c Delta))/(8 (Delta^2 + 1)) - ((theta^2_c - 3 theta_c + 3) e^(theta_c))/(16) dot \
dot Re[E_1(theta_c - i theta_c Delta)] + ((theta_c^2 + 3 theta_c + 3) e^(-theta_c))/(16) Re[E_1(-theta_c - i 
theta_c Delta)]} +\
+ (4 omega^2_"max")/(E^2_0 lambda) lr(\{, size: #300%) -E^2_0/omega^3_"max" (pi/2 e^(-theta_c) - 1/2[e^(-theta_c) 
Re[E_1(theta_c - i theta_c Delta)] + e^(theta_c) Re[E_1(-theta_c - i theta_c Delta)]]) - \
- E_0/omega^2_"max" ((Delta sin(theta_c Delta))/(2 (Delta^2 + 1)) + (theta_c - 1)/4 e^(theta_c) 
Re[E_1(theta_c - i theta_c Delta)] + (theta_c + 1)/4 e^(-theta_c) dot \
dot Re[E_1(-theta_c - i theta_c Delta)]) - 1/omega_"max" ((Delta (3 Delta^2 + 5) sin(theta_c Delta))/(8 
(Delta^2 + 1)^2) + \
+ (theta_c cos(theta_c Delta))/(8 (Delta^2 + 1)) - ((theta^2_c - 3 theta_c + 3)e^(theta_c))/16 dot
dot Re[E_1(theta_c - i theta_c Delta)] + \
+ ((theta^2_c + 3 theta_c + 3) e^(-theta_c))/(16) 
Re[E_1(- theta_c - i theta_c Delta)]) + E^3_0/omega^3_"max" dot \ 
dot (pi/(2 c) e^(-theta_c c) - 1/(2 c) [e^(-theta_c c) Re[E_1(theta_c c - i theta_c Delta)] + e^(theta_c c) 
Re[E_1(-theta_c c - i theta_c Delta)]])lr(}, size: #300%) - \
- (4 omega^3_"max")/(E^2_0 lambda^2) {-1/omega_"max" ((Delta sin(theta_c Delta))/(2(Delta^2 + 1)) + 
(theta_c cos(theta_c Delta))/(8(Delta^2 + 1)) - ((theta_c^2 - 3theta_c + 3)e^(theta_c))/(16) dot \
dot Re[E_1(theta_c - i theta_c Delta)] + ((theta^2_c - 3theta_c + 3)e^(-theta_c))/(16) 
Re[E_1(-theta_c - i theta_c Delta)]) + sqrt(E_0)/(omega^(3/2)_"max") dot \
dot (pi/(2 a) e^(-theta_c a) -1/(2a)[e^(-theta_c a) Re[E_1(theta_c a - i theta_c Delta)]
+ e^(theta_c a) Re[E_1(-theta_c a) - i theta_c Delta]] + \
+ pi/(2 b) e^(-theta_c b) - 1/(2 b) [e^(-theta_c b) Re[E_1(theta_c b - i theta_c Delta)] + e^(theta_c b) 
Re[E_1(- theta_c b - i theta_c Delta)]])}}
$

