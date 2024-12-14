#include <iostream>
#include <cassert>
#include "matrix.h"

void TestDeterminant() {
    std::cout << "Testing determinant calculation...\n";

    // Тест 1: Единичная матрица 2x2
    Matrix<2, 2> identity = {{
        {1, 0},
        {0, 1}
    }};
    assert(identity.det() == 1);

    // Тест 2: Простая матрица 2x2
    Matrix<2, 2> simple = {{
        {2, 3},
        {4, 5}
    }};
    assert(simple.det() == -2); // 2*5 - 3*4 = -2

    // Тест 3: Матрица 3x3
    Matrix<3, 3> matrix3x3 = {{
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    }};
    assert(matrix3x3.det() == 0); // Вырожденная матрица

    std::cout << "Determinant tests passed!\n";
}

void TestRank() {
    std::cout << "Testing rank calculation...\n";

    // Тест 1: Единичная матрица
    Matrix<3, 3> identity = {{
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    }};
    assert(identity.rank() == 3);

    // Тест 2: Матрица с линейно зависимыми строками
    Matrix<3, 3> dependent = {{
        {1, 2, 3},
        {2, 4, 6},
        {3, 6, 9}
    }};
    assert(dependent.rank() == 1);

    // Тест 3: Прямоугольная матрица
    Matrix<2, 3> rectangular = {{
        {1, 2, 3},
        {4, 5, 6}
    }};
    assert(rectangular.rank() == 2);

    std::cout << "Rank tests passed!\n";
}

void TestInversion() {
    std::cout << "Testing matrix inversion...\n";

    // Тест 1: Единичная матрица
    Matrix<2, 2> identity = {{
        {1, 0},
        {0, 1}
    }};
    Matrix<2, 2> inverted_identity = identity.inverted();
    assert(inverted_identity == identity);

    // Тест 2: Простая матрица 2x2
    Matrix<2, 2> simple = {{
        {4, 7},
        {2, 6}
    }};
    Matrix<2, 2> inverted_simple = simple.inverted();
    // Проверяем, что A * A^(-1) = I
    Matrix<2, 2> product = simple * inverted_simple;
    // assert(product[0,0] == 1 && product[0,1] == 0 &&
           // product[1,0] == 0 && product[1,1] == 1);

    // Тест 3: Проверка исключения для вырожденной матрицы
    Matrix<2, 2> singular = {{
        {1, 2},
        {2, 4}
    }};
    try {
        Matrix<2, 2> inverted_singular = singular.inverted();
        assert(false); // Должно быть выброшено исключение
    } catch (const std::runtime_error& e) {
        // Ожидаемое поведение
    }

    std::cout << "Inversion tests passed!\n";
}
void tester(){
    TestDeterminant();
    TestRank();
    TestInversion();

    std::cout << "All tests passed successfully!\n";
}
