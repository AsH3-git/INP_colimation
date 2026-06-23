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
$
J_2 = limits(integral)^Delta_(-infinity) (Im(e^(i theta_c psi)) psi)/(psi^2 + c^2) d psi = 
Im limits(integral)^Delta_(-infinity) (e^(i theta_c psi) psi)/(psi^2 + c^2) d psi = K \
K = limits(integral)^infinity_(-infinity) - limits(integral)^infinity_Delta = "In"_1 -"In"_2 \
"In"_1 = limits(integral)^infinity_(-infinity) (e^(i theta_c psi) psi)/ (psi^2 + c^2) d psi \
"Контур - верхняя полуокружность." z_0 = i c - "полюс первого порядка" \
limits(Res(f(psi)))_(psi = z_0) = phi(z_0)/(zeta'(z_0)) = (e^(i theta_c psi) psi)/(2)|_(z=z_0) = e^(- theta_c c)/2 \
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
