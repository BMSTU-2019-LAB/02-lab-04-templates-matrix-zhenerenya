// Copyright 2020  zhenerenya zhenyachmchk@gmail.com

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_
#include <math.h>
#include <stdio.h>

#include <limits>
#include <type_traits>

template <class T>
class Matrix {
  int rows;     //строки
  int columns;  //столбцы
  T** m;        //матрица типа T

 public:
  //констуктор пустого объекта
  Matrix();
  //конструктор создания объекта
  Matrix(int rows, int columns);
  //конструктор копирования марицы
  Matrix(const Matrix<T>& s);
  ~Matrix();
  //инты
  int get_rows() const;
  int get_columns() const;
  //лонг
  size_t Rows() const;
  static bool is_comparible(const Matrix<T>& s1, const Matrix<T>& s2) {
    return ((s1.rows == s2.rows) && (s1.columns == s2.columns));
  }
  //перегрузка операции присваивание
  Matrix<T>& operator=(const Matrix<T>& s);

  //операция сложения матриц
  Matrix<T> operator+(const Matrix<T>& s2);

  //операция вычитание матриц
  Matrix<T> operator-(const Matrix<T>& s2);

  //перегрузка операции ображения по индексу
  T* operator[](int index) const;

  //умножение матриц
  Matrix<T> operator*(const Matrix<T>& s2);

  //опредилитель
  T determinant(const Matrix<T>& s);

  //вычеркивание срок
  Matrix<T> deleted(const Matrix<T>& s, int x, int y);

  //матрица алгебраических дополнений
  Matrix<T> minors(const Matrix<T>& s);

  //транспонирование
  Matrix<T> transp(const Matrix<T>& s);

  //умножение матрицы на число
  void mult(const T& alfa, Matrix<T>& s);

  //нахождение обратной матрицы
  Matrix<T> Inverse();

  //сравнение матриц
  template <class v>
  friend bool operator==(const Matrix<v>& s1, const Matrix<v>& s2);

  template <class v>
  friend bool operator!=(const Matrix<v>& s1, const Matrix<v>& s2);
};

//конструктор без параметров
template <class T>
Matrix<T>::Matrix() {
  rows = 0;
  columns = 0;
  m = nullptr;
}

//конструктор создания объекта
template <class T>
Matrix<T>::Matrix(int rows, int columns) {
  this->rows = rows;
  this->columns = columns;
  m = new T*[rows];
  for (int i = 0; i < rows; i++) {
    m[i] = new T[columns];
  }
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      m[i][j] = 0;
    }
  }
}

//конструктор копирования марицы
template <class T>
Matrix<T>::Matrix(const Matrix<T>& s) {
  this->rows = s.rows;
  this->columns = s.columns;
  m = new T*[rows];
  for (int i = 0; i < rows; i++) {
    m[i] = new T[columns];
  }
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; ++j) {
      if (s[i][j])
        m[i][j] = s[i][j];
      else
        m[i][j] = 0;
    }
  }
}

template <class T>
Matrix<T>::~Matrix() {
  if (m != nullptr) {
    for (int i = 0; i < rows; i++) {
      delete[] m[i];
    }
    delete[] m;
  }
}

//инты
template <class T>
int Matrix<T>::get_rows() const {
  return rows;
}
template <class T>
int Matrix<T>::get_columns() const {
  return columns;
}

//лонг
template <class T>
size_t Matrix<T>::Rows() const {
  return rows;
}

//перегрузка операции присваивание
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& s) {
  for (int i = 0; i < this->rows; ++i) {
    delete[] m[i];
  }
  delete[] m;
  this->rows = s.rows;
  this->columns = s.columns;
  m = new T*[rows];
  for (int i = 0; i < rows; i++) {
    m[i] = new T[columns];
  }
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; ++j) {
      m[i][j] = s.m[i][j];
    }
  }
  return *this;
}

//операция сложения матриц
template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& s2) {
  Matrix<T> s1(rows, columns);  //сумма таблиц
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      s1[i][j] = m[i][j];
    }
  }
  if (is_comparible(s1, s2)) {
    Matrix<T> new_matrix(rows, columns);  // можно ли сложить
    new_matrix = Matrix(s1.get_rows(),
                        s1.get_columns());  //создали матрицу нужного размера
    for (int i = 0; i < s1.get_rows(); ++i) {
      for (int j = 0; j < s1.get_columns(); ++j) {
        new_matrix[i][j] = s1[i][j] + s2[i][j];
      }
    }
    return new_matrix;
  } else {
    return {};
  }
}

//операция вычитание матриц
template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& s2) {
  Matrix<T> s1(rows, columns);  //сумма таблиц
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      s1[i][j] = m[i][j];
    }
  }
  if (is_comparible(s1, s2)) {
    Matrix<T> new_matrix(rows, columns);  // можно ли сложить
    new_matrix = Matrix(s1.get_rows(),
                        s1.get_columns());  //создали матрицу нужного размера
    for (int i = 0; i < s1.get_rows(); ++i) {
      for (int j = 0; j < s1.get_columns(); ++j) {
        new_matrix[i][j] = s1[i][j] - s2[i][j];
      }
    }
    return new_matrix;
  } else {
    return {};
  }
}

//перегрузка операции ображения по индексу
template <class T>
T* Matrix<T>::operator[](int index) const {  //не изменяет значения
  return m[index];
}

//умножение матриц
template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& s2) {
  Matrix<T> s1(rows, columns);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < rows; ++j) {
      s1[i][j] = m[i][j];
    }
  }
  if (s1.columns == s2.rows) {
    Matrix<T> new_matrix(rows, columns);
    new_matrix = Matrix(s1.rows, s2.columns);
    int index_row = 0;
    while (index_row < s1.rows) {
      int index_column = 0;
      while (index_column < s2.columns) {
        T sum = 0;
        for (int k = 0; k < s1.columns; ++k) {
          sum += s1[index_row][k] * s2[k][index_column];
        }
        new_matrix[index_row][index_column] = sum;
        index_column++;
      }
      index_row++;
    }
    return new_matrix;
  } else
    return {};
}

//опредилитель
template <class T>
T Matrix<T>::determinant(const Matrix<T>& s) {
  T det = 0;
  if (s.columns == 1)
    det = s[0][0];
  else {
    if (s.columns == 2)
      det = s[0][0] * s[1][1] - s[0][1] * s[1][0];
    else {

        for (int j = 0; j < s.columns; ++j) {
          if ((j % 2)==1)
            det +=s[0][j]* (-1) * determinant(deleted(s, 0, j));
          else
            det += s[0][j]*determinant(deleted(s, 0, j));
        }

    }
  }
  return det;
}

//вычеркивание срок
template <class T>
Matrix<T> Matrix<T>::deleted(const Matrix<T>& s, int x, int y) {
  Matrix<T> new_matrix(s.rows - 1, s.columns - 1);
  int new_rows = 0;
  int new_columns = 0;
  for (int i = 0; i < s.rows; i++) {
    if (i != x) {
      for (int j = 0; j < s.rows; ++j) {
        if (j != y) {
          new_matrix[new_rows][new_columns] = s[i][j];
          ++new_columns;
        }
      }
      ++new_rows;
      new_columns = 0;
    }
  }
  return new_matrix;
}

//матрица алгебраических дополнений
template <class T>
Matrix<T> Matrix<T>::minors(const Matrix<T>& s) {
  Matrix<T> matrix_minors;
  matrix_minors = Matrix(s.rows, s.columns);
  for (int i = 0; i < s.rows; i++) {
    for (int j = 0; j < s.columns; j++) {
      if ((i + j) % 2 ==1)
        matrix_minors[i][j] = (-1) * determinant(deleted(s, i, j));
      else
        matrix_minors[i][j] = determinant(deleted(s, i, j));
    }
  }
  return matrix_minors;
}

//транспонирование
template <class T>
Matrix<T> Matrix<T>::transp(const Matrix<T>& s) {
  Matrix<T> newm(rows, columns);
  for (int i = 0; i < rows ; i++) {
    for (int j = 0; j < rows ; j++) {
      newm[j][i] = s[i][j];
      newm[i][j] = s[j][i];
    }
  }
  return newm;
}

//умножение матрицы на число
template <class T>
void Matrix<T>::mult(const T& alfa, Matrix<T>& s) {
  for (int i = 0; i < s.rows; ++i) {
    for (int j = 0; j < s.columns; ++j) {
      s[i][j] *= alfa;
    }
  }
}

//нахождение обратной матрицы
template <class T>
Matrix<T> Matrix<T>::Inverse() {
  //определитель x
  //нахождение минора x
  //ставим его на место x
  //транспонируем матрицу x
  //умножаем на 1/опредиоитель
  if (this->columns == this->rows) {
    Matrix<T> new_matrix(columns, columns);
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < rows; ++j) {
        new_matrix[i][j] = m[i][j];
      }
    }
    int det = determinant(new_matrix);
    new_matrix = transp(minors(new_matrix));
     mult(1. / det, new_matrix);
    return new_matrix;
  } else
    return {};
}

//сравнение матриц
template <class T>
bool operator==(const Matrix<T>& s1, const Matrix<T>& s2) {
  if (s1.get_rows() == s2.get_rows() && s1.get_columns() == s2.get_columns()) {
    if (std::is_floating_point<T>::value) {
      for (int i = 0; i < s1.rows; ++i) {
        for (int j = 0; j < s1.columns; ++j) {
          if (fabs(s1[i][j] - s2[i][j]) > std::numeric_limits<double>::epsilon())
            return false;
        }
      }
    }

    for (int i = 0; i < s1.rows; ++i) {
      for (int j = 0; j < s1.columns; ++j) {
        if (s1[i][j] != s2[i][j]) return false;
      }
    }
    return true;
  } else
    return false;
}

template <class T>
bool operator!=(const Matrix<T>& s1, const Matrix<T>& s2) {
  if (s1 == s2)
    return false;
  else
    return true;
}

#endif  // INCLUDE_MATRIX_HPP_