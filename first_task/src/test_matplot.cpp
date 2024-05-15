#include "matplot/matplot.h"

int main() {
  // matplot::backend("Agg");
  // Создаем векторы для x и y 
  std::vector<double> x = matplot::linspace(0, 10, 100);
  std::vector<double> y = x;

  // Строим график 
  matplot::plot(x, y);

  // Настройка осей и заголовка 
  matplot::xlabel("x");
  matplot::ylabel("y");
  matplot::title("График y=x");

  // Отображение графика
  matplot::save("/home/alex/parproga/first_task/src/graph.png"); 
  return 0;
}