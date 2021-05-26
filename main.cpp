// Gauspos.cpp : ���� ���� �������� ������� "main". ����� ���������� � ������������� ���������� ���������.
//
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>

#include <omp.h>

int* pSerialPivotPos; //  ������ �����, �������� ��  ���������
int* pSerialPivotIter; // ��������, �� ���� ���� ���� ��������
// ������� ��� ������ ����������� ������� �� ��������� ��������
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
// ������� ��������� ����������� ������� �� ��������� ��������
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
// ������� �������� ���'�� �� ���������� �������� ��'����
void ProcessInitialization(double*& pMatrix, double*& pVector,
	double*& pResult, int& Size) {
	// ������������ ������ ������� �� �������
	do {
		printf("\nEnter size of the matrix and the vector: ");
		scanf("%d", &Size);
		printf("\nChosen size = %d \n", Size);
		if (Size <= 0)
			printf("\nSize of objects must be greater than 0!\n");
	} while (Size <= 0);
	// �������� ���'��
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];
	// ����������� ������� �� ��������� ��������
	DummyDataInitialization(pMatrix, pVector, Size);
	//RandomDataInitialization(pMatrix, pVector, Size);
}
// ������� ��� ������ �������
void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	int i, j;
	for (i = 0; i < RowCount; i++) {
		for (j = 0; j < ColCount; j++)
			printf("%7.4f ", pMatrix[i * RowCount + j]);
		printf("\n");
	}
}
// ������� ��� ������ ��������
void PrintVector(double* pVector, int Size) {
	int i;
	for (i = 0; i < Size; i++)
		printf("%7.4f ", pVector[i]);
}

// ����� �������� �����
int FindPivotRow(double* pMatrix, int Size, int Iter) {
	int PivotRow = -1;
	// ������ �������� �����
	int MaxValue = 0; // ��������  �������� �������� 
	int i;

	// �������� �����, � ����� ���������� ������������ �������
	for (i = 0; i < Size; i++) {
		if ((pSerialPivotIter[i] == -1) &&
			(fabs(pMatrix[i * Size + Iter]) > MaxValue)) {
			PivotRow = i;
			MaxValue = fabs(pMatrix[i * Size + Iter]);
		}
	}
	return PivotRow;
}

// ��������� �������
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

// �������� ��������
void SerialGaussianElimination(double* pMatrix, double* pVector, int Size) {
	int Iter;

	int PivotRow;
	// ����� ��������� �������� �����
	for (Iter = 0; Iter < Size; Iter++) {
		// ����� �������� �����
		PivotRow = FindPivotRow(pMatrix, Size, Iter);
		pSerialPivotPos[Iter] = PivotRow;

		pSerialPivotIter[PivotRow] = Iter;
		SerialColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
	}
}
// ��������� ���
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
// ������� ��� ��������� ��������� �����
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
// ������� ���������� ��������������� �������
void ProcessTermination(double* pMatrix, double* pVector, double* pResult)
{
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
}
void main() {
	double* pMatrix; // ������� ����� �������
	double* pVector;	// ���� ������� ����� �������
	double* pResult; // ������ ����������
	int Size;	// ������ ��������� ������� �� �������
	time_t start, finish;
	double duration;
	printf("Serial Gauss algorithm for solving linear systems\n");
	// �������� ���'�� �� ���������� �������� ��'����
	ProcessInitialization(pMatrix, pVector, pResult, Size);


	//printf("Initial Matrix \n");
//PrintMatrix(pMatrix, Size, Size);
	printf("Initial Vector \n");
	PrintVector(pVector, Size);
	// ��������� ��������� �����
	start = clock();
	SerialResultCalculation(pMatrix, pVector, pResult, Size);
	finish = clock();
	duration = (finish - start) / double(CLOCKS_PER_SEC);

	// ���� ������� ����������
	printf("\n Result Vector: \n");
	PrintVector(pResult, Size);

	printf("\n Time of execution: %f\n", duration);

	// ���������� ��������������� �������
	ProcessTermination(pMatrix, pVector, pResult);

}
// ������ ���������: CTRL+F5 ��� ���� "�������" > "������ ��� �������"
// ������� ���������: F5 ��� ���� "�������" > "��������� �������"

// ������ �� ������ ������ 
//   1. � ���� ������������ ������� ����� ��������� ����� � ��������� ���.
//   2. � ���� Team Explorer ����� ������������ � ������� ���������� ��������.
//   3. � ���� "�������� ������" ����� ������������� �������� ������ ������ � ������ ���������.
//   4. � ���� "������ ������" ����� ������������� ������.
//   5. ��������������� �������� ������ ���� "������" > "�������� ����� �������", ����� ������� ����� ����, ��� "������" > "�������� ������������ �������", ����� �������� � ������ ������������ ����� ����.
//   6. ����� ����� ������� ���� ������ �����, �������� ������ ���� "����" > "�������" > "������" � �������� SLN-����.
