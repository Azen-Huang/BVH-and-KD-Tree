#include<iostream>
#include<stdlib.h>
#include<math.h>

class matrix3
{
public:
	matrix3(float a11, float a12, float a13,
		   float a21, float a22, float a23,
		   float a31, float a32, float a33);
	matrix3();
	~matrix3();

	//矩陣相乘
	static matrix3 mutiply(matrix3 A,matrix3 B);
	//矩陣相加
	static matrix3 add(matrix3 A,matrix3 B);
	//矩陣的逆
	static matrix3 inverse(matrix3 A);

	//二維陣列，儲存矩陣元素
	float M[3][3];
};

matrix3::matrix3()
{
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
			this->M[i][j] = 0.00;
	}
}

matrix3::matrix3(
			   float a11, float a12, float a13,
			   float a21, float a22, float a23,
			   float a31, float a32, float a33)
{
	M[0][0] = a11; M[0][1] = a12; M[0][2] = a13;
	M[1][0] = a21; M[1][1] = a22; M[1][2] = a23; 
	M[2][0] = a31; M[2][1] = a32; M[2][2] = a33; 
}

matrix3::~matrix3()
{

}

///Function:兩個矩陣相乘
///para:A-矩陣A
///para:B-矩陣B
///return：矩陣AB
matrix3 matrix3::mutiply(matrix3 A, matrix3 B)
{
	matrix3 Result;
	for(int m=0;m<3;m++)
	{  
		for(int s=0;s<3;s++)
		{
			for(int n=0;n<3;n++)
			{  
				Result.M[m][s]+=A.M[m][n]*B.M[n][s];  
			}  
		}  
	}
	return Result;
}

///Function:兩個矩陣相加
///para:A-矩陣A
///para:B-矩陣B
///return：矩陣AB
matrix3 matrix3::add(matrix3 A, matrix3 B)
{
	matrix3 Result;
	for(int m=0;m<3;m++)
	{  
		for(int s=0;s<3;s++)
		{
			Result.M[m][s]=A.M[m][s]+B.M[m][s]; 
		}  
	}
	return Result;
}

///Function:求矩陣的逆
///para:A-矩陣A
///return：矩陣A的逆矩陣
matrix3 matrix3::inverse(matrix3 A)
{
	float detA = 
		(A.M[0][0] * A.M[1][1] * A.M[2][2]) + 
		(A.M[1][0] * A.M[2][1] * A.M[0][2]) +
		(A.M[0][1] * A.M[1][2] * A.M[2][0]) -
		(A.M[0][2] * A.M[1][1] * A.M[2][0]) -
		(A.M[1][2] * A.M[2][1] * A.M[0][0]) -
		(A.M[1][0] * A.M[0][1] * A.M[2][2]);
	/*
	00 01 02
	10 11 12
	20 21 22
	*/
	matrix3 Result;
	Result.M[0][0] = 1 / detA * (A.M[1][1] * A.M[2][2] - A.M[2][1] * A.M[1][2]);
	Result.M[1][0] = -1 / detA * (A.M[1][0] * A.M[2][2] - A.M[2][0] * A.M[1][2]);
	Result.M[2][0] = 1 / detA * (A.M[1][0] * A.M[2][1] - A.M[2][0] * A.M[1][1]);

	Result.M[0][1] = -1 / detA * (A.M[0][1] * A.M[2][2] - A.M[2][1] * A.M[0][2]);
	Result.M[1][1] = 1 / detA * (A.M[0][0] * A.M[2][2] - A.M[2][0] * A.M[0][2]);
	Result.M[2][1] = -1 / detA * (A.M[0][0] * A.M[2][1] - A.M[2][0] * A.M[0][1]);

	Result.M[0][2] = 1 / detA * (A.M[0][1] * A.M[1][2] - A.M[1][1] * A.M[0][2]);
	Result.M[1][2] = -1 / detA * (A.M[0][0] * A.M[1][2] - A.M[1][0] * A.M[0][2]);
	Result.M[2][2] = 1 / detA * (A.M[0][0] * A.M[1][1] - A.M[1][0] * A.M[0][1]);

	return Result;
}