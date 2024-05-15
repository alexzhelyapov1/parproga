#include <matplot/matplot.h>

int main() {
  using namespace matplot;

  // Создаем данные для графика
  std::vector<double> x = linspace(0, 10, 100);
  std::vector<double> y = x;

  // Строим график
  plot(x, y);

  // Задаем название осей
  xlabel("x");
  ylabel("y");

  // Сохраняем график в файл "plot.png"
  save("plot.png");

  return 0;
}