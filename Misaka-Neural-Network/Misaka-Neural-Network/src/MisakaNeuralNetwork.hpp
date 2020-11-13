#ifndef MISAKA_NEURAL_NETWORK_H
#define MISAKA_NEURAL_NETWORK_H

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

namespace Misaka
{

    class Matrix
    {
    public:
        int mRows, mCols;

    private:
        std::vector<std::vector<double>> mData;

    public:
        /**
         * @brief       Construct a new Matrix object.
        */
        Matrix() : mRows(0), mCols(0) { /* Empty */ }

        /**
         * @brief           Construct a new Matrix object
         * @param pRows     Number of mRows in matrix.
         * @param pCols     Number of columns in matrix.
         * @return          void
        */
        Matrix(int pRows, int pCols) : mRows(pRows), mCols(pCols)
        {
            for (int i = 0; i < mRows; i++)
            {
                std::vector<double> row;
                for (int j = 0; j < mCols; j++)
                {
                    row.push_back(0);
                }
                mData.push_back(row);
            }
        }

        /**
         * @brief           Construct a new Matrix object
         * @param pData     Matrix data.
         * @example         Creating 2x3 matrix:
         *                  Matrix mat = Matrix({{ 1, 2, 3 }, { 4, 5, 6 }});
         * @return          void
        */
        Matrix(const std::vector<std::vector<double>>& pData) : mRows(pData.size()), mCols(pData[0].size())
        {
            mData = pData;
        }

        /**
         * @brief           Multiplies this matrix with scalar.
         * @param pValue    Scalar.
         * @return          void
        */
        void Mult(double pValue)
        {
            for (auto& row : mData)
            {
                for (auto& value : row)
                {
                    value *= pValue;
                }
            }
        }

        /**
         * @brief           Multiplies this matrix with other matrix.
         * @param pMatrix   Other matrix object.
         * @return          void
        */
        void Mult(Matrix& pMatrix)
        {
            // Dot procut.
            if (mCols == pMatrix.mRows)
            {
                Matrix resultMatrix(mRows, pMatrix.mCols);
                for (int i = 0; i < resultMatrix.mRows; i++)
                {
                    for (int j = 0; j < resultMatrix.mCols; j++)
                    {
                        double sum = 0;
                        for (int k = 0; k < mCols; k++)
                        {
                            sum += mData[i][k] * pMatrix[k][j];
                        }
                        resultMatrix[i][j] = sum;
                    }
                }
                mData = resultMatrix.mData;
                mRows = resultMatrix.mRows;
                mCols = resultMatrix.mCols;
            }
            // Element wise product.
            else if (mCols == pMatrix.mCols && mRows == pMatrix.mRows)
            {
                for (int i = 0; i < mRows; i++)
                {
                    for (int j = 0; j < mCols; j++)
                    {
                        mData[i][j] = mData[i][j] * pMatrix[i][j];
                    }
                }
            }
            // Error.
            else
            {
                std::cout << "Matrix multiplication error!" << std::endl;
            }
        }

        /**
         * @brief           Multiplies this matrix with scalar.
         * @param pValue    Scalar.
         * @return          Matrix
        */
        Matrix operator *(double pValue)
        {
            Matrix matrix = *this;
            matrix.Mult(pValue);
            return matrix;
        }

        /**
         * @brief           Multiplies this matrix with other matrix.
         * @param pMatrix   Other matrix object.
         * @return          Matrix
        */
        Matrix operator *(Matrix& pMatrix)
        {
            Matrix matrix = *this;
            matrix.Mult(pMatrix);
            return matrix;
        }

        /**
         * @brief           Add value to every element in matrix.
         * @param pValue    Value to add.
         * @return          void
        */
        void Add(double pValue)
        {
            for (auto& row : mData)
            {
                for (auto& value : row)
                {
                    value += pValue;
                }
            }
        }

        /**
         * @brief           Element wise matrix addition.
         * @param pMatrix   Other matrix object.
         * @return          void
        */
        void Add(Matrix& pMatrix)
        {
            for (int i = 0; i < mRows; i++)
            {
                for (int j = 0; j < mCols; j++)
                {
                    mData[i][j] = mData[i][j] + pMatrix[i][j];
                }
            }
        }

        /**
         * @brief           Add value to every element in matrix.
         * @param pValue    Value to add.
         * @return          Matrix
        */
        Matrix operator +(double pValue)
        {
            Matrix matrix = *this;
            matrix.Add(pValue);
            return matrix;
        }

        /**
         * @brief           Element wise matrix addition.
         * @param pMatrix   Other matrix object.
         * @return          Matrix
        */
        Matrix operator +(Matrix& pMatrix)
        {
            Matrix matrix = *this;
            matrix.Add(pMatrix);
            return matrix;
        }

        /**
         * @brief           Subtract value from every element in matrix.
         * @param pValue    Value to subtract.
         * @return          void
        */
        void Sub(double pValue)
        {
            for (auto& row : mData)
            {
                for (auto& value : row)
                {
                    value -= pValue;
                }
            }
        }

        /**
         * @brief           Element wise matrix subtraction.
         * @param pMatrix   Other matrix object.
         * @return          void
        */
        void Sub(Matrix& pMatrix)
        {
            for (int i = 0; i < mRows; i++)
            {
                for (int j = 0; j < mCols; j++)
                {
                    mData[i][j] = mData[i][j] - pMatrix[i][j];
                }
            }
        }

        /**
         * @brief           Subtract value from every element in matrix.
         * @param pValue    Value to subtract.
         * @return          Matrix
        */
        Matrix operator -(double pValue)
        {
            Matrix matrix = *this;
            matrix.Sub(pValue);
            return matrix;
        }

        /**
         * @brief           Element wise matrix subtraction.
         * @param pMatrix   Other matrix object.
         * @return          Matrix
        */
        Matrix operator -(Matrix& pMatrix)
        {
            Matrix matrix = *this;
            matrix.Sub(pMatrix);
            return matrix;
        }

        /**
         * @brief           Set random values for matrix elements.
         * @param min       Mininum value.
         * @param max       Maximum value.
        */
        void Randomize(double min = -1, double max = 1)
        {
            for (auto& row : mData)
            {
                for (auto& value : row)
                {
                    value = GetRandomNumber(min, max);
                }
            }
        }

        /**
         * @brief           Transpose this matrix object.
         * @return          void
        */
        void Transpose()
        {
            Matrix matrix(mCols, mRows);
            for (int i = 0; i < mRows; i++)
            {
                for (int j = 0; j < mCols; j++)
                {
                    matrix[j][i] = mData[i][j];
                }
            }

            mData = matrix.mData;
            mRows = matrix.mRows;
            mCols = matrix.mCols;
        }

        /**
         * @brief           Creates matrix equal to provided transpose matrix. Does not change provided matrix.
         * @param pMatrix   Other matrix object.
         * @return          Matrix
        */
        static Matrix Transpose(Matrix& pMatrix)
        {
            Matrix matrix(pMatrix.mCols, pMatrix.mRows);
            for (int i = 0; i < pMatrix.mRows; i++)
            {
                for (int j = 0; j < pMatrix.mCols; j++)
                {
                    matrix[j][i] = pMatrix.mData[i][j];
                }
            }
            return matrix;
        }

        /**
         * @brief           Returns matrix row.
         * @param i         Row number.
         * @return          std::vector<dobule>&
        */
        std::vector<double>& operator [](int i)
        {
            return mData[i];
        }

        /**
         * @brief           Print matrix values to standard output.
         * @return          void
        */
        void Print()
        {
            int num = 0;
            std::cout << "[ Matrix Data ]\n";
            for (const auto& row : mData)
            {
                std::cout << "\t" << num++ << "\t:\t";
                for (const auto& value : row)
                {
                    std::cout << std::fixed << value << "\t";
                }
                std::cout << "\n";
            }
        }

    private:

        /**
         * @brief           Generate and return random number in range.
         * @param min       Minimum value.
         * @param max       Maximum value.
         * @return          double
        */
        double GetRandomNumber(double min, double max)
        {
            return min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (max - min)));
        }
    };


    class NeuralNetwork
    {
    private:
        std::map<std::string, Matrix> mWeights;
        std::map<std::string, Matrix> mBiases;
        double mLearningRate = 0.1;

    public:
        /**
         * @brief           Construct a new NeuralNetwork object.
         * @param pInput    Number of input nodes.
         * @param pHidden   Number of hidden nodes.
         * @param pOutput   Number of output nodes.
        */
        NeuralNetwork(int pInput, int pHidden, int pOutput)
        {
            (mWeights["IH"] = Matrix(pHidden, pInput)).Randomize(-1, 1);
            (mWeights["HO"] = Matrix(pOutput, pHidden)).Randomize(-1, 1);

            (mBiases["H"] = Matrix(pHidden, 1)).Randomize(-1, 1);
            (mBiases["O"] = Matrix(pOutput, 1)).Randomize(-1, 1);
        }

        void SetLearningRate(double pLearningRate)
        {
            mLearningRate = pLearningRate;
        }

        /**
         * @brief           Generate output from provided input data.
         * @param pData     Input data.
         * @return          Matrix
        */
        Matrix Query(std::vector<double> pData)
        {
            Matrix input({ pData });
            input.Transpose();

            Matrix hidden = mWeights["IH"] * input;
            hidden.Add(mBiases["H"]);
            ApplySigmoid(hidden);

            Matrix output = mWeights["HO"] * hidden;
            output.Add(mBiases["O"]);
            ApplySigmoid(output);

            return output;
        }

        /**
         * @brief           Train neural network.
         * @param pInput    Input data.
         * @param pTarget   Target data.
         * @return          void
        */
        void Train(std::vector<double> pInput, std::vector<double> pTarget)
        {
            Matrix input({ pInput });
            input.Transpose();

            Matrix hidden = mWeights["IH"] * input;
            hidden.Add(mBiases["H"]);
            ApplySigmoid(hidden);

            Matrix output = mWeights["HO"] * hidden;
            output.Add(mBiases["O"]);
            ApplySigmoid(output);

            Matrix target = Matrix({ pTarget });
            target.Transpose();

            Matrix outputErrors = target - output;
            Matrix outputGradients = ApplySigmoidDerivative(output);

            outputGradients.Mult(outputErrors);
            outputGradients.Mult(mLearningRate);

            Matrix hiddenT = Matrix::Transpose(hidden);
            Matrix whod = outputGradients * hiddenT;
            mWeights["HO"].Add(whod);
            mBiases["O"].Add(outputGradients);

            Matrix whot = Matrix::Transpose(mWeights["HO"]);

            Matrix hiddenErrors = whot * outputErrors;
            Matrix hiddenGradients = ApplySigmoidDerivative(hidden);

            hiddenGradients.Mult(hiddenErrors);
            hiddenGradients.Mult(mLearningRate);

            Matrix inputT = Matrix::Transpose(input);
            Matrix wihd = hiddenGradients * inputT;
            mWeights["IH"].Add(wihd);
            mBiases["H"].Add(hiddenGradients);
        }

    private:
        /**
         * @brief           Apply sigmoid function to every element in matrix.
         * @param pMatrix   Matrix object.
         * @return void
        */
        void ApplySigmoid(Matrix& pMatrix)
        {
            for (int i = 0; i < pMatrix.mRows; i++)
            {
                for (int j = 0; j < pMatrix.mCols; j++)
                {
                    pMatrix[i][j] = 1.0 / (1.0 + std::exp(-pMatrix[i][j]));
                }
            }
        }

        /**
         * @brief           Apply sigmoid derivative to every element in matrix.
         * @param pMatrix   Matrix object.
         * @return          Matrix
        */
        Matrix ApplySigmoidDerivative(Matrix& pMatrix)
        {
            Matrix matrix = Matrix(pMatrix.mRows, pMatrix.mCols);
            for (int i = 0; i < pMatrix.mRows; i++)
            {
                for (int j = 0; j < pMatrix.mCols; j++)
                {
                    matrix[i][j] = pMatrix[i][j] * (1.0 - pMatrix[i][j]);
                }
            }
            return matrix;
        }
    };
}

#endif