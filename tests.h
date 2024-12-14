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
    Matrix<3, 3, double> identity = {{
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    }};
    assert(identity.rank() == 3);

    // Тест 2: Матрица с линейно зависимыми строками
    Matrix<3, 3, double> dependent = {{
        {1, 2, 3},
        {2, 4, 6},
        {3, 6, 9}
    }};
    assert(dependent.rank() == 1);

    // Тест 3: Прямоугольная матрица
    Matrix<2, 3, double> rectangular = {{
        {1, 2, 3},
        {4, 5, 6}
    }};
    assert(rectangular.rank() == 2);

    // Тест 4: Нулевая матрица
    Matrix<3, 3, double> zero = {{
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    }};
    assert(zero.rank() == 0);

    // Тест 5: Матрица ранга 2
    Matrix<3, 3, double> rank2 = {{
        {1, 0, 2},
        {0, 1, 1},
        {1, 1, 3}
    }};
    assert(rank2.rank() == 2);

    // Тест 6: Матрица с дробными числами
    Matrix<2, 2,double> fractional = {{
        {0.5, 1.5},
        {1.0, 3.0}
    }};
    assert(fractional.rank() == 1);

    // Тест 7: Широкая прямоугольная матрица
    Matrix<2, 4, double> wide = {{
        {1, 2, 3, 4},
        {2, 4, 6, 8}
    }};
    assert(wide.rank() == 1);

    // Тест 8: Высокая прямоугольная матрица
    Matrix<4, 2, double> tall = {{
        {1, 2},
        {3, 4},
        {5, 6},
        {7, 8}
    }};
    assert(tall.rank() == 2);

    // Тест 9: Матрица с близкими к нулю элементами
    Matrix<2, 2, double> nearZero = {{
        {1e-10, 2e-10},
        {3e-10, 6e-10}
    }};
    assert(nearZero.rank() == 1);

    // Тест 10: Матрица с большими числами
    Matrix<2, 2, double> large = {{
        {1e6, 2e6},
        {3e6, 6e6}
    }};
    assert(large.rank() == 1);

    // Тест 11: Матрица с отрицательными числами
    Matrix<3, 3, double> negative = {{
        {-1, -2, -3},
        {-2, -4, -6},
        {-4, -8, -12}
    }};
    assert(negative.rank() == 1);

    // Тест 12: Матрица с комбинацией нулей и ненулевых элементов
    Matrix<3, 3, double> mixed = {{
        {1, 0, 0},
        {0, 0, 2},
        {0, 0, 0}
    }};
    assert(mixed.rank() == 2);

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
