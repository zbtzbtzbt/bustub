//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// p0_starter.h
//
// Identification: src/include/primer/p0_starter.h
//
// Copyright (c) 2015-2020, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include "common/exception.h"
#include "common/logger.h"
namespace bustub {
    using std::vector;

/**
 * The Matrix type defines a common
 * interface for matrix operations.
 */
template <typename T>
class Matrix {
protected:
    /**
     * TODO(P0): Add implementation [Done]
     *
     * Construct a new Matrix instance.
     * @param rows The number of rows
     * @param cols The number of columns
     *
     */
    Matrix(int rows, int cols) {
        rows_ = rows;
        cols_ = cols;
        for (int i = 0; i < rows_; i++) {
            vector<T> tmp_row(cols_);
            matrix_.push_back(tmp_row);
        }
        linear_ =  new T[rows_*cols_+1];
    }
    
    /** The number of rows in the matrix */
    int rows_;
    /** The number of columns in the matrix */
    int cols_;
    /** Matrix */
    vector <vector<T>> matrix_;
    T* linear_;
    /**
     * TODO(P0): Allocate the array in the constructor.
     * TODO(P0): Deallocate the array in the destructor.
     * A flattened array containing the elements of the matrix.
     */
    
 public:
    // =0 after function means pure virtual function
    /** @return The number of rows in the matrix */
    virtual int GetRowCount() const = 0;
    
    /** @return The number of columns in the matrix */
    virtual int GetColumnCount() const = 0;
    
    /**
     * Get the (i,j)th matrix element.
     *
     * Throw OUT_OF_RANGE if either index is out of range.
     *
     * @param i The row index
     * @param j The column index
     * @return The (i,j)th matrix element
     * @throws OUT_OF_RANGE if either index is out of range
     */
    virtual T GetElement(int i, int j) const = 0;
    
    /**
     * Set the (i,j)th matrix element.
     *
     * Throw OUT_OF_RANGE if either index is out of range.
     *
     * @param i The row index
     * @param j The column index
     * @param val The value to insert
     * @throws OUT_OF_RANGE if either index is out of range
     */
    virtual void SetElement(int i, int j, T val) = 0;
    
    /**
     * Fill the elements of the matrix from `source`.
     *
     * Throw OUT_OF_RANGE in the event that `source`
     * does not contain the required number of elements.
     *
     * @param source The source container
     * @throws OUT_OF_RANGE if `source` is incorrect size
     */
    virtual void FillFrom(const std::vector <T>& source) = 0;
    
    /**
     * Destroy a matrix instance.
     * TODO(P0): Add implementation
     */
    virtual ~Matrix() = default;
};

/**
 * The RowMatrix type is a concrete matrix implementation.
 * It implements the interface defined by the Matrix type.
 */
template <typename T>
class RowMatrix : public Matrix<T> {
public:
    /**
     * TODO(P0): Add implementation [Done]
     *
     * Construct a new RowMatrix instance.
     * @param rows The number of rows
     * @param cols The number of columns
     */
//    explicit RowMatrix(int rows, int cols) : Matrix<T>(rows,cols){
//        // use base class to init
//        //LOG_INFO("RowMatrix(%d,%d) use base class to init",rows,cols);
//        data_ = &linear_;
//    }
    RowMatrix(int rows,int cols): Matrix<T>(rows, cols){
        data_ = &(Matrix<T>::linear_);
    }
    
    
    /**
     * TODO(P0): Add implementation [Done]
     * @return The number of rows in the matrix
     */
    int GetRowCount() const override { return Matrix<T>::rows_; }
    
    /**
     * TODO(P0): Add implementation [Done]
     * @return The number of columns in the matrix
     */
    int GetColumnCount() const override { return Matrix<T>::cols_; }
    
    /**
     * TODO(P0): Add implementation [Done]
     *
     * Get the (i,j)th matrix element.
     *
     * Throw OUT_OF_RANGE if either index is out of range.
     *
     * @param i The row index
     * @param j The column index
     * @return The (i,j)th matrix element
     * @throws OUT_OF_RANGE if either index is out of range
     */
    T GetElement(int i, int j) const override {
        if (i < 0 || i > Matrix<T>::rows_ || j < 0 || j > Matrix<T>::cols_) {
            //throw NotImplementedException{"RowMatrix::GetElement(i,j) not implemented."};
            //LOG_INFO("GetElement(%d, %d) out of range!",i,j);
            T* tmp = new T;
            return *tmp;
        }
        return Matrix<T>::matrix_[i][j];
    }
    
    /**
     * Set the (i,j)th matrix element.
     *
     * Throw OUT_OF_RANGE if either index is out of range.
     *
     * @param i The row index
     * @param j The column index
     * @param val The value to insert
     * @throws OUT_OF_RANGE if either index is out of range
     */
    void SetElement(int i, int j, T val) override {
        if (i < 0 || i > Matrix<T>::rows_ || j < 0 || j > Matrix<T>::cols_) {
            //throw NotImplementedException{"RowMatrix::GetElement(i,j) not implemented."};
            //LOG_INFO("SetElement(%d, %d) out of range!",i,j);
            return ;
        }
        Matrix<T>::matrix_[i][j] = val;
        return;
    }
    
    /**
     * TODO(P0): Add implementation
     *
     * Fill the elements of the matrix from `source`.
     *
     * Throw OUT_OF_RANGE in the event that `source`
     * does not contain the required number of elements.
     *
     * @param source The source container
     * @throws OUT_OF_RANGE if `source` is incorrect size
     */
    void FillFrom(const std::vector <T>& source) override {
        throw NotImplementedException{"RowMatrix::FillFrom() not implemented."};
    }
    
    /**
     * TODO(P0): Add implementation
     *
     * Destroy a RowMatrix instance.
     */
    ~RowMatrix() override = default;

private:
    /**
     * A 2D array containing the elements of the matrix in row-major format.
     *
     * TODO(P0):
     * - Allocate the array of row pointers in the constructor.
     * - Use these pointers to point to corresponding elements of the `linear` array.
     * - Don't forget to deallocate the array in the destructor.
     */
    T** data_;
};

/**
 * The RowMatrixOperations class defines operations
 * that may be performed on instances of `RowMatrix`.
 */
template <typename T>
class RowMatrixOperations {
public:
    /**
     * Compute (`matrixA` + `matrixB`) and return the result.
     * Return `nullptr` if dimensions mismatch for input matrices.
     * @param matrixA Input matrix
     * @param matrixB Input matrix
     * @return The result of matrix addition
     */
    static std::unique_ptr <RowMatrix<T>> Add(const RowMatrix<T>* matrixA, const RowMatrix<T>* matrixB) {
        // TODO(P0): Add implementation
        
        if (matrixA->GetColumnCount() != matrixB->GetColumnCount() ||
            matrixA->GetRowCount() != matrixB->GetRowCount()) {
            //LOG_INFO("Add size mismatch",i,j);
            return std::unique_ptr<RowMatrix<T>>(nullptr);
        }
        const int rows = matrixA->GetRowCount();
        const int cols = matrixA->GetColumnCount();
    
        RowMatrix<T>* matrixC = new RowMatrix<T>(rows,cols);
    
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrixC->SetElement(i, j,matrixA->GetElement(i, j) + matrixB->GetElement(i, j));
            }
        }
        
        return std::unique_ptr<RowMatrix<T>>(matrixC);
    }
    
    /**
     * Compute the matrix multiplication (`matrixA` * `matrixB` and return the result.
     * Return `nullptr` if dimensions mismatch for input matrices.
     * @param matrixA Input matrix
     * @param matrixB Input matrix
     * @return The result of matrix multiplication
     */
    static std::unique_ptr <RowMatrix<T>> Multiply(const RowMatrix<T>* matrixA, const RowMatrix<T>* matrixB) {
        // TODO(P0): Add implementation
        const int c_rows_ = matrixA->GetColumnCount();
        const int c_cols = matrixB->GetRowCount();
        
        RowMatrix<T>* matrixC = new RowMatrix<T>(c_rows_,c_cols);
        return std::unique_ptr<RowMatrix<T>>(matrixC);
    }
    
    /**
     * Simplified General Matrix Multiply operation. Compute (`matrixA` * `matrixB` + `matrixC`).
     * Return `nullptr` if dimensions mismatch for input matrices.
     * @param matrixA Input matrix
     * @param matrixB Input matrix
     * @param matrixC Input matrix
     * @return The result of general matrix multiply
     */
    static std::unique_ptr <RowMatrix<T>> GEMM(const RowMatrix<T>* matrixA, const RowMatrix<T>* matrixB,
                                               const RowMatrix<T>* matrixC) {
        // TODO(P0): Add implementation
        RowMatrix<T>* matrixD =Add(Multiply(matrixA,matrixB),matrixC);
        return std::unique_ptr<RowMatrix<T>>(matrixD);
    }
};
}  // namespace bustub
