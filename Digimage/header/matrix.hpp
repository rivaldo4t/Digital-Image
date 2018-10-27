#pragma once
#include <iostream>
#include <vector>

class Vector
{
private:
	std::vector<float> vec;
public:
	Vector(std::vector<float> v = std::vector<float>(3)) : vec(v) {}
	Vector(const float x, const float y, const float z) 
	{ 
		vec = std::vector<float>(3);
		vec[0] = x; vec[1] = y; vec[2] = z; 
	}
	Vector(const Vector& v) 
	{ 
		vec = v.getVector(); 
	}
	
	const std::vector<float>& getVector() const
	{ return vec; }

	Vector& operator = (const Vector& v)
	{
		vec = v.getVector();
		return *this;
	}

	float& operator [] (const int i) { return vec[i]; }
	const float& operator [] (const int i) const { return vec[i]; }

	friend std::ostream& operator << (std::ostream &os, const Vector& v)
	{ return os << v[0] << " " << v[1] << " " << v[2] << std::endl; }

	void normalize()
	{
		vec[0] /= vec[2];
		vec[1] /= vec[2];
	}
};

class Matrix
{
protected:
	std::vector<std::vector<float>> mat;
public:
	Matrix(std::vector<std::vector<float>> m = std::vector<std::vector<float>>(3, std::vector<float>(3))) : mat(m) {}
	Matrix(const float m00, const float m01, const float m02, const float m10, const float m11, 
		const float m12, const float m20, const float m21, const float m22)
	{
		mat = std::vector<std::vector<float>>(3, std::vector<float>(3));
		mat[0][0] = m00; mat[0][1] = m01; mat[0][2] = m02; mat[1][0] = m10; mat[1][1] = m11; 
		mat[1][2] = m12; mat[2][0] = m20; mat[2][1] = m21; mat[2][2] = m22;
	}
	Matrix(const Matrix& m)
	{
		mat = m.getMatrix(); 
	}

	const std::vector<std::vector<float>>& getMatrix() const
	{ return mat; }

	Matrix& operator = (const Matrix& m)
	{
		mat = m.getMatrix();
		return *this;
	}

	const float& operator () (const int i, const int j) const { return mat[i][j]; }

	friend std::ostream& operator << (std::ostream &os, const Matrix& m)
	{
		return os	<< m(0, 0) << " " << m(0, 1) << " " << m(0, 2) << std::endl
					<< m(1, 0) << " " << m(1, 1) << " " << m(1, 2) << std::endl
					<< m(2, 0) << " " << m(2, 1) << " " << m(2, 2) << std::endl;
	}

	const Vector operator * (const Vector& v) const
	{
		std::vector<float> mult;
		for (int i = 0; i < 3; ++i)
		{
			float sum = 0;
			for (int j = 0; j < 3; ++j)
			{
				sum += mat[i][j] * v[j];
			}
			mult.push_back(sum);
		}
		return Vector(mult);
	}

	const Matrix operator * (const Matrix& m) const
	{
		std::vector<std::vector<float>> multMat;
		for (int i = 0; i < 3; ++i)
		{
			std::vector<float> multVec;
			for (int j = 0; j < 3; ++j)
			{
				float sum = 0;
				for (int k = 0; k < 3; ++k)
				{
					sum += mat[i][k] * m(k, j);
				}
				multVec.push_back(sum);
			}
			multMat.push_back(multVec);
		}
		return Matrix(multMat);
	}

	virtual void invertMatrix() {}
};

class RotationMatrix : public Matrix
{
private:
	float theta;
public:
	RotationMatrix(float t) : Matrix(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), theta(t) {}
	void invertMatrix()
	{
		mat[0][1] = -sin(theta);
		mat[1][0] = sin(theta);
	}
};

class ScaleMatrix : public Matrix
{
private:
	float sx, sy;
public:
	ScaleMatrix(float x, float y) : Matrix(x, 0, 0, 0, y, 0, 0, 0, 1), sx(x), sy(y) {}
	void invertMatrix()
	{
		mat[0][0] = 1 / sx;
		mat[1][1] = 1 / sy;
	}
};

class TranslateMatrix : public Matrix
{
private:
	float tx, ty;
public:
	TranslateMatrix(float x, float y) : Matrix(1, 0, tx, 0, 1, ty, 0, 0, 1), tx(x), ty(y) {}
	void invertMatrix()
	{
		mat[0][2] = -tx;
		mat[1][2] = -ty;
	}
};

class ShearMatrix : public Matrix
{
private:
	float shx, shy;
public:
	ShearMatrix(float x, float y) : Matrix(1, x, 0, y, 1, 0, 0, 0, 1), shx(x), shy(y) {}
	void invertMatrix()
	{
		mat[0][1] = -shx;
		mat[1][0] = -shy;
		mat[2][2] = 1 - shx*shy;
	}
};

class PerspectiveMatrix : public Matrix
{
private:
	float px, py;
public:
	PerspectiveMatrix(float x, float y) : Matrix(1, 0, 0, 0, 1, 0, x, y, 1), px(x), py(y) {}
	void invertMatrix()
	{
		mat[2][0] = -px;
		mat[2][1] = -py;
	}
};