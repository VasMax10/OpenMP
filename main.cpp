// Gauspos.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>

#include <omp.h>

int* pSerialPivotPos; //  вудучі рядків, вибраних на  ітераціях
int* pSerialPivotIter; // Ітерації, на яких ряди були ведучими
// Функція для простої ініціалізації матриці та векторних елементів
void DummyDataInitialization(double* pMatrix, double* pVector, int Size) {
	int i, j;
	for (i = 0; i < Size; i++) {
		pVector[i] = i + 1;
		for (j = 0; j < Size; j++) {
			if (j <= i)
				pMatrix[i * Size + j] = 1;
			else
				pMatrix[i * Size + j] = 0;
		}
	}
}
// Функція випадкової ініціалізації матриці та векторних елементів
void RandomDataInitialization(double* pMatrix, double* pVector, int Size)
{
	int i, j;
	srand(unsigned(clock()));
	for (i = 0; i < Size; i++) {
		pVector[i] = rand() / double(1000);
		for (j = 0; j < Size; j++) {
			if (j <= i)
				pMatrix[i * Size + j] = rand() / double(1000);
			else
				pMatrix[i * Size + j] = 0;
		}
	}
}
// Функція виділення пам'яті та визначення елементів об'єктів
void ProcessInitialization(double*& pMatrix, double*& pVector,
	double*& pResult, int& Size) {
	// Встановлення розміру матриці та вектора
	do {
		printf("\nEnter size of the matrix and the vector: ");
		scanf("%d", &Size);
		printf("\nChosen size = %d \n", Size);
		if (Size <= 0)
			printf("\nSize of objects must be greater than 0!\n");
	} while (Size <= 0);
	// Виділення пам'яті
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];
	// Ініціалізація матриці та векторних елементів
	DummyDataInitialization(pMatrix, pVector, Size);
	//RandomDataInitialization(pMatrix, pVector, Size);
}
// Функція для виводу матриці
void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	int i, j;
	for (i = 0; i < RowCount; i++) {
		for (j = 0; j < ColCount; j++)
			printf("%7.4f ", pMatrix[i * RowCount + j]);
		printf("\n");
	}
}
// Функція для виводу векторна
void PrintVector(double* pVector, int Size) {
	int i;
	for (i = 0; i < Size; i++)
		printf("%7.4f ", pVector[i]);
}

// Пошук ведучого рядка
int FindPivotRow(double* pMatrix, int Size, int Iter) {
	int PivotRow = -1;
	// Індекс ведучого рядка
	int MaxValue = 0; // Значення  ведучого елемента 
	int i;

	// Вибераємо рядок, в якому зберігається максимальний елемент
	for (i = 0; i < Size; i++) {
		if ((pSerialPivotIter[i] == -1) &&
			(fabs(pMatrix[i * Size + Iter]) > MaxValue)) {
			PivotRow = i;
			MaxValue = fabs(pMatrix[i * Size + Iter]);
		}
	}
	return PivotRow;
}

// занулення стовпця
void SerialColumnElimination(double* pMatrix, double* pVector, int Pivot,
	int Iter, int Size) {
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[Pivot * Size + Iter];
	#pragma omp parallel for default(none) shared(PivotValue, PivotFactor)
	for (int i = 0; i < Size; i++) {
		if (pSerialPivotIter[i] == -1) {
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue;
			for (int j = Iter; j < Size; j++) {
				pMatrix[i * Size + j] -= PivotFactor * pMatrix[Pivot * Size + j];
			}
			pVector[i] -= PivotFactor * pVector[Pivot];
		}
	}
}

// Гауссова елімінація
void SerialGaussianElimination(double* pMatrix, double* pVector, int Size) {
	int Iter;

	int PivotRow;
	// Номер поточного ведучого рядка
	for (Iter = 0; Iter < Size; Iter++) {
		// Пошук ведучого рядка
		PivotRow = FindPivotRow(pMatrix, Size, Iter);
		pSerialPivotPos[Iter] = PivotRow;

		pSerialPivotIter[PivotRow] = Iter;
		SerialColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
	}
}
// обернений хід
void SerialBackSubstitution(double* pMatrix, double* pVector,
	double* pResult, int Size) {
	int RowIndex, Row;
	#pragma omp parallel for default(none) shared(RowIndex, Row)
	for (int i = Size - 1; i >= 0; i--) {
		RowIndex = pSerialPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
		for (int j = 0; j < i; j++) {
			Row = pSerialPivotPos[j];
			pVector[j] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}
// Функція для виконання алгоритму Гауса
void SerialResultCalculation(double* pMatrix, double* pVector,
	double* pResult, int Size) {

	pSerialPivotPos = new int[Size];
	pSerialPivotIter = new int[Size];
	for (int i = 0; i < Size; i++) {
		pSerialPivotIter[i] = -1;
	}

	SerialGaussianElimination(pMatrix, pVector, Size);

	SerialBackSubstitution(pMatrix, pVector, pResult, Size);

	delete[] pSerialPivotPos;
	delete[] pSerialPivotIter;
}
// Функція завершення обчислювального процесу
void ProcessTermination(double* pMatrix, double* pVector, double* pResult)
{
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
}
void main() {
	double* pMatrix; // Матриця лінійної системи
	double* pVector;	// Праві частини лінійної системи
	double* pResult; // Вектор результату
	int Size;	// Розміри початкової матриці та вектора
	time_t start, finish;
	double duration;
	printf("Serial Gauss algorithm for solving linear systems\n");
	// Виділення пам'яті та визначення елементів об'єктів
	ProcessInitialization(pMatrix, pVector, pResult, Size);


	//printf("Initial Matrix \n");
//PrintMatrix(pMatrix, Size, Size);
	printf("Initial Vector \n");
	PrintVector(pVector, Size);
	// Виконання алгоритму Гауса
	start = clock();
	SerialResultCalculation(pMatrix, pVector, pResult, Size);
	finish = clock();
	duration = (finish - start) / double(CLOCKS_PER_SEC);

	// Друк вектора результату
	printf("\n Result Vector: \n");
	PrintVector(pResult, Size);

	printf("\n Time of execution: %f\n", duration);

	// Припинення обчислювального процесу
	ProcessTermination(pMatrix, pVector, pResult);

}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
