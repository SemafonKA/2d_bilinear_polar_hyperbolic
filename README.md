# МКЭ для гиперболического уравнения в полярной системе координат, схема трёхслойная явная, базисные функции билинейные на прямоугольниках.

[toc]

## О курсовой

Данная курсовая работа тесно связана с [прошлой курсовой работой](https://github.com/SemafonKA/2D_bilinear_polar_ellipse), и по факту всё, что меняется - добавляется сетка по времени.

Решаемое уравнение в общем виде: 

$$
-{\rm div}(\lambda {\rm grad}(u)) + \gamma u + \sigma \frac{\partial u}{\partial t} + \chi \frac{\partial^{2}u}{\partial t^{2}} = f
$$

Область интегрирования: $\Omega$

Граница интегрирования: $S = S_1 \lor S_2 \lor S_3$

Начальные условия:

$$
\begin{gather*}
   \left. u \right|_{t=t_{0}} = u_{0},  \\
   \left. \frac{\partial u}{\partial t} \right|_{t=t_{0}} = u^{'}_{0}, 
\end{gather*}
$$

и тогда следующий слой $u_{1}$ будет считаться по следующей формуле:

$$
u_{1} \approx u_{0} + \left. \frac{\partial u}{\partial t} \right|_{t=t_{0}} \cdot (t-t_{0}),
$$

либо можем задать следующий слой сразу:

$$
\begin{gather*}
   \left. u \right|_{t=t_{0}} = u_{0}, \\
   \left. u \right|_{t=t_{1}} = u_{1}.
\end{gather*}
$$

Краевые условия:

$$
\begin{gather*}
   u|_{s_1} = u_{g} \\
   \lambda \frac{\partial u}{\partial n} = \theta \\
   \lambda \frac{\partial u}{\partial n} |_{s_3} + \beta(u|_{s_3} - u_{\beta}) = 0
\end{gather*}
$$

$\lambda$ - коэффициент диффузии,

$\beta$ - коэффициент теплообмена.

Согласно условиям варианта, в задаче используется полярная система координат $(r, \varphi)$. Для неё формулы операторов уравнения определяются следующим видом:

$$
{\rm grad}\space v = (\frac{\partial v}{\partial r}, \frac{1}{r} \frac{\partial v}{\partial \varphi})
$$

$$
{\rm div}\space \overrightarrow{Q} = \frac{1}{r} \frac{\partial (rQ_1)}{\partial r} + \frac{1}{r} \frac{\partial Q_2}{\partial \varphi}
$$

А само уравнение в таком случае будет выглядеть следующим образом:

$$
-\frac{1}{r}\frac{\partial}{\partial r}\left(r\lambda \frac{\partial u}{\partial r}\right) - \frac{1}{r^{2}}\frac{\partial}{\partial \varphi}\left(\lambda \frac{\partial u}{\partial \varphi}\right) + \gamma u + \sigma \frac{\partial u}{\partial t} + \chi \frac{\partial^{2}u}{\partial t^{2}} = f
$$

В учебнике конечномерная аппроксимация строится для случая, когда $\gamma = 0$, и формула выглядит так:

$$ 
-{\rm div}(\lambda {\rm grad}(u)) + \sigma \frac{\partial u}{\partial t} + \chi \frac{\partial^{2}u}{\partial t^{2}} = f
$$

И такая формула с учётом системы координат примет следующий вид:

$$
-\frac{1}{r}\frac{\partial}{\partial r}\left(r\lambda \frac{\partial u}{\partial r}\right) - \frac{1}{r^{2}}\frac{\partial}{\partial \varphi}\left(\lambda \frac{\partial u}{\partial \varphi}\right) + \sigma \frac{\partial u}{\partial t} + \chi \frac{\partial^{2}u}{\partial t^{2}} = f
$$

Подробности + тесты находятся в файле [Курсач-6-семак.md](Курсач-6-семак.md), а также в [отчёте](<Самсонов%20С.%20ПМ-01%20курсовая.pdf>).

## Тестирование

В целом, тестирование не сильно отличается от прошлой курсовой работы. Только теперь нас практически не интересуют тесты по пространству, а больше интересуют тесты по временной сетке. Примеры тестов находятся [здесь](Курсач-6-семак.md#тестирование).

Не стоит делать тесты при маленьком $\chi$, поскольку тогда задача будет сводиться к эллиптической, и почему-то плохо решаться. Крч не стоит так делать.

## Сборка

Для сборки потребуется установленный cmake версии &ge; 3.0.0, а также компилятор c++ с поддержкой стандарта языка c++20. 

### Если используется компилятор Visual Studio:

1. Запустить консоль разработчика Visual Studio (developer powershell for VS 20XX);
2. Перейти в папку с CMakeLists.txt файлом (`cd "\*путь до директории\*"`);
3. Перейти в папку build, либо если её нет - создать и перейти (`cd build`, `mkdir build`);
4. Сгенерировать решение при помощи команды `cmake ..`;
5. Запустить проект Visual Studio и сделать сборку через неё, либо прописать в консоль `cmake --build .`

### Если используется компилятор gcc или clang:

К сожалению, особо не шарю за них, но скорее всего сборка будет идентичной, за исключением запуска консоли разработчика (это касается только компилятора Visual Studio).

