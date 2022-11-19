#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() { this->init_matrix(3, 3); }

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::invalid_argument("Invalid matrix size");
  }

  this->init_matrix(rows, cols);
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  this->init_matrix(other._rows, other._cols);
  this->copy_matrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  this->init_matrix(other._rows, other._cols);
  this->copy_matrix(other);
  other.remove_matrix();
}

S21Matrix::~S21Matrix() { this->remove_matrix(); }

bool S21Matrix::eq_matrix(const S21Matrix& other) {
  bool result = FAILURE;
  if (this->_cols == other._cols && this->_rows == other._rows) {
    bool flag = 1;
    for (int i = 0; i < this->_rows; i++) {
      for (int j = 0; j < this->_cols; j++) {
        if (fabs(this->_matrix[i][j] - other._matrix[i][j]) > 0.000001) {
          result = FAILURE;
          flag = 0;
        } else {
          if (flag) {
            result = SUCCESS;
          }
        }
      }
    }
  }
  return result;
}

void S21Matrix::sum_matrix(const S21Matrix& other) {
  if (this->_rows != other._rows || this->_cols != other._cols) {
    throw std::invalid_argument("Matrix sizes are not equal");
  }

  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_matrix[i][j] += other._matrix[i][j];
    }
  }
}

void S21Matrix::sub_matrix(const S21Matrix& other) {
  if (this->_rows != other._rows || this->_cols != other._cols) {
    throw std::invalid_argument("Matrix sizes are not equal");
  }

  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_matrix[i][j] -= other._matrix[i][j];
    }
  }
}

void S21Matrix::mul_number(const double num) {
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_matrix[i][j] *= num;
    }
  }
}

void S21Matrix::mul_matrix(const S21Matrix& other) {
  if (this->_cols != other._rows) {
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix");
  }

  S21Matrix result(this->_rows, other._cols);

  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      for (int k = 0; k < other._rows; k++) {
        result._matrix[i][j] += this->_matrix[i][k] * other._matrix[k][j];
      }
    }
  }
  *this = result;
}

S21Matrix S21Matrix::transpose() {
  S21Matrix result(this->_cols, this->_rows);

  for (int i = 0; i < this->_cols; i++) {
    for (int j = 0; j < this->_rows; j++) {
      result._matrix[i][j] = this->_matrix[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::calc_complements() {
  if (this->_cols != this->_rows || this->_rows < 2) {
    throw std::invalid_argument("The matrix is not square");
  }

  S21Matrix result(this->_rows, this->_cols);

  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      S21Matrix minor_matrix = this->matrix_minor(i, j);
      result._matrix[i][j] = minor_matrix.determinant() * pow(-1, i + j);
    }
  }
  return result;
}

double S21Matrix::determinant() {
  if (this->_cols != this->_rows) {
    throw std::invalid_argument("The matrix is not square");
  }

  double determinant = 0;
  if (this->_rows == 1) {
    determinant = this->_matrix[0][0];
  } else if (this->_rows == 2) {
    determinant = this->_matrix[0][0] * this->_matrix[1][1] -
                  this->_matrix[1][0] * this->_matrix[0][1];
  } else {
    for (int k = 0; k < this->_rows; k++) {
      S21Matrix minor(this->_rows - 1, this->_cols - 1);
      for (int i = 1; i < this->_rows; i++) {
        for (int j = 0; j < this->_rows; j++) {
          if (j == k)
            continue;
          else if (j < k)
            minor._matrix[i - 1][j] = this->_matrix[i][j];
          else
            minor._matrix[i - 1][j - 1] = this->_matrix[i][j];
        }
      }
      determinant += pow(-1, k + 2) * this->_matrix[0][k] * minor.determinant();
    }
  }
  return determinant;
}

S21Matrix S21Matrix::inverse_matrix() {
  double det = this->determinant();
  if (det == 0 || this->_cols != this->_rows) {
    throw std::invalid_argument(
        "The matrix is not square or the determinant is zero");
  }

  S21Matrix alg = this->calc_complements();
  S21Matrix trans = alg.transpose();
  trans.mul_number(1 / det);
  return trans;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.sum_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.sub_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.mul_matrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.mul_number(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) {
  return this->eq_matrix(other);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    this->remove_matrix();
    this->init_matrix(other._rows, other._cols);
    this->copy_matrix(other);
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->sum_matrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->sub_matrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->mul_matrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  this->mul_number(num);
  return *this;
}

double& S21Matrix::operator()(int row, int col) const {
  if (row >= this->_rows || row < 0 || col >= this->_cols || col < 0) {
    throw std::invalid_argument("Indices outside the matrix");
  }

  return this->_matrix[row][col];
}

int S21Matrix::GetRows() const { return this->_rows; }

int S21Matrix::GetCols() const { return this->_cols; }

void S21Matrix::SetRows(int num) {
  if (num <= 0) {
    throw std::invalid_argument("Impossible meaning");
  }

  S21Matrix new_matrix(num, this->_cols);
  int count_rows = 0;

  if (this->_rows <= num) {
    count_rows = this->_rows;
  } else {
    count_rows = num;
  }

  for (int i = 0; i < count_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      new_matrix._matrix[i][j] = this->_matrix[i][j];
    }
  }

  *this = new_matrix;
}

void S21Matrix::SetCols(int num) {
  if (num <= 0) {
    throw std::invalid_argument("Impossible meaning");
  }

  S21Matrix new_matrix(this->_rows, num);
  int count_cols = 0;

  if (this->_cols <= num) {
    count_cols = this->_cols;
  } else {
    count_cols = num;
  }

  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < count_cols; j++) {
      new_matrix._matrix[i][j] = this->_matrix[i][j];
    }
  }

  *this = new_matrix;
}

S21Matrix S21Matrix::matrix_minor(int rows, int cols) {
  S21Matrix result(this->_rows - 1, this->_cols - 1);

  int ci = 0, cj = 0;
  for (int i = 0; i < this->_rows; i++) {
    if (i == rows) continue;
    for (int j = 0; j < this->_cols; j++) {
      if (j == cols) continue;
      result._matrix[ci][cj] = this->_matrix[i][j];
      cj++;
    }
    ci++;
    cj = 0;
  }
  return result;
}

void S21Matrix::copy_matrix(const S21Matrix& other) {
  this->_rows = other._rows;
  this->_cols = other._cols;
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_matrix[i][j] = other._matrix[i][j];
    }
  }
}

void S21Matrix::remove_matrix() {
  if (this->_cols != -1 || this->_rows != -1) {
    for (int i = 0; i < this->_rows; i++) {
      delete[] this->_matrix[i];
    }
    delete[] this->_matrix;
    this->_rows = -1;
    this->_cols = -1;
  }
}

void S21Matrix::init_matrix(int rows, int cols) {
  this->_rows = rows;
  this->_cols = cols;
  this->_matrix = new double*[this->_rows];
  for (int i = 0; i < this->_rows; i++) {
    this->_matrix[i] = new double[this->_cols];
  }
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_matrix[i][j] = 0;
    }
  }
}
