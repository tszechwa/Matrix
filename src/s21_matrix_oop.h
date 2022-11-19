#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#define SUCCESS 1
#define FAILURE 0

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  int _rows, _cols;
  double** _matrix;

 public:
  S21Matrix();  // Конструктор по умолчанию +!
  S21Matrix(int rows,
            int cols);  // Конструктор с количеством строк и столбцов +!
  S21Matrix(const S21Matrix& other);  // Конструктор копирования +!
  S21Matrix(S21Matrix&& other);  // Конструктор переноса +!
  ~S21Matrix();                  // Деструктор +!

  bool eq_matrix(const S21Matrix& other);  // Равенство матриц +!
  void sum_matrix(const S21Matrix& other);  // Сложение матриц(искл.:различная
                                            // размерность матриц) +!
  void sub_matrix(const S21Matrix& other);  // Разность матриц(искл.:различная
                                            // размерность матриц) +!
  void mul_number(const double num);  // Умножение матрицы на число +!
  void mul_matrix(
      const S21Matrix& other);  // Умножение двух матриц(искл.:число столбцов 1
                                // матрицы != числу строк 2 матрицы) +!
  S21Matrix transpose();  // Транспонирование матрицы +!
  S21Matrix calc_complements();  // Матрица алгебраический дополнений(искл.:не
                                 // является квадратной) +!
  double determinant();          // Определитель матрицы +!
  S21Matrix inverse_matrix();  // Обратная матрица(искл.:определитель = 0) +!

  S21Matrix operator+(
      const S21Matrix&
          other);  // Сложение матриц(искл.:различная размерность матриц) +!
  S21Matrix operator-(
      const S21Matrix&
          other);  // Разность матриц(искл.:различная размерность матриц) +!
  S21Matrix operator*(
      const S21Matrix& other);  // Умножение двух матриц(искл.:число столбцов 1
                                // матрицы != числу строк 2 матрицы) +!
  S21Matrix operator*(const double num);  // Умножение матрицы на число +!
  bool operator==(const S21Matrix& other);  // Проверка на равенство +!
  S21Matrix& operator=(const S21Matrix& other);  // Присваивание +!
  S21Matrix& operator+=(
      const S21Matrix&
          other);  // Присвоение сложения(искл.:различная размерность матриц) +!
  S21Matrix& operator-=(
      const S21Matrix&
          other);  // Присвоение разности(искл.:различная размерность матриц) +!
  S21Matrix& operator*=(
      const S21Matrix&
          other);  // Присвоение умножения двух матриц(искл.:число столбцов 1
                   // матрицы != числу строк 2 матрицы) +!
  S21Matrix& operator*=(
      const double num);  // Присвоение умножения матрицы на число +!
  double& operator()(int row, int col)
      const;  // Индексация(искл.:индекс за пределами матрицы) +!

  int GetRows() const;    // +!
  int GetCols() const;    // +!
  void SetRows(int num);  // +!
  void SetCols(int num);  // +!

  void copy_matrix(const S21Matrix& other);
  void remove_matrix();
  void init_matrix(int rows, int cols);
  S21Matrix matrix_minor(int rows, int cols);
};

#endif  // SRC_S21_MATRIX_OOP_H_
