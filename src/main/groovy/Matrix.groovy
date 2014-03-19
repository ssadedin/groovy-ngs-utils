import groovy.transform.CompileStatic;

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.RealMatrix;

class MatrixColumn implements Iterable {
    
    int columnIndex
    
    Matrix sourceMatrix
    
    Object getAt(Object index) {
        if(index instanceof Integer)
            sourceMatrix.dataRef[index][columnIndex]
        else
        if(index instanceof List)
            sourceMatrix.dataRef[index].collect { it[columnIndex] }
    }
    
    int size() {
        sourceMatrix.matrix.columnDimension
    }
    
    
    Iterator iterator() {
        return new Iterator<Double>() {
            int index = 0
            int max = MatrixColumn.this.size()
            
            boolean hasNext() {
                return index<=max
            }
            
            Double next() {
                return sourceMatrix.dataRef[index++][columnIndex]
            }
            
            void remove() { throw new UnsupportedOperationException() }
        }
    }
    
    String toString() {
        "[" + this.collect {it}.join(",") + "]"
    }
}

class Matrix {
    
    @Delegate
    Array2DRowRealMatrix matrix

    public Matrix(int rows, int columns) {
        matrix = new Array2DRowRealMatrix(rows, columns)
    }
    
    @CompileStatic
    public Matrix(int rows, int columns, List<Double> data) {
        matrix = new Array2DRowRealMatrix(rows, columns)
        int i=0
        for(int r=0; r<rows; ++r) {
            for(int c=0; c<columns;++c) {
                matrix.dataRef[r][c] = (double)data[i++]
            }
        }
    }
    
    @CompileStatic
    public Matrix(int rows, int columns, double[] matrixData) {
        matrix = new Array2DRowRealMatrix(rows, columns)
        int i=0
        for(int r=0; r<rows; ++r) {
            for(int c=0; c<columns;++c) {
                matrix.dataRef[r][c] = matrixData[++i]
            }
        }
    }
      
    public Matrix(Array2DRowRealMatrix m) {
        matrix = m
    }
     
    @CompileStatic
    MatrixColumn col(int n) {
//        matrix.getColumn(n)
        new MatrixColumn(columnIndex:n, sourceMatrix: this)
    }
    
    List<MatrixColumn> getColumns() {
        (0..<matrix.columnDimension).collect { col(it) }
    }
    
    @CompileStatic
    double [] row(int n) {
        matrix.getRow(n)
    }
    
    @CompileStatic
    Object getAt(Object n) {
        if(n == null) {
            return getColumns()
        }
        else
        if(n instanceof Integer)
            return matrix.dataRef[(int)n]
        else
        if(n instanceof List) {
            List l = (List)n
            if(l.size() == 0) // Seems to happen with m[][2] type syntax
                return getColumns()
            else
                return matrix.dataRef[(int)l[0]][(int)l[1]]
        }
    }
    
    @CompileStatic
    def collect(Closure c) {
        List<Object> results = new ArrayList(matrix.dataRef.size())
        int rowIndex = 0;
        if(c.maximumNumberOfParameters == 1) {
            for(double [] row in matrix.dataRef) {
                results.add(c(row))
            }
        }
        else 
        if(c.maximumNumberOfParameters == 2) {
            for(double [] row in matrix.dataRef) {
                results.add(c(row, rowIndex))
                ++rowIndex
            }
        }
        return results
    }    
    
    
    @CompileStatic
    Matrix grep(Closure c) {
        List<Integer> keepRows = []
        int rowIndex = 0;
        if(c.maximumNumberOfParameters == 1) {
            for(double [] row in matrix.dataRef) {
                if(c(row))
                    keepRows.add(rowIndex)
                ++rowIndex
            }
        }
        else 
        if(c.maximumNumberOfParameters == 2) {
            for(double [] row in matrix.dataRef) {
                if(c(row, rowIndex))
                    keepRows.add(rowIndex)
                ++rowIndex
            }
        }
        
        return new Matrix(new Array2DRowRealMatrix(matrix.data[keepRows]))
    }    
    
    /**
     * Transforms a matrix by processing each element through the given
     * closure. The closure must take either one argument or three arguments.
     * The one argument version is only passed data values, while the 
     * three argument version is passed the data value and also the row and column
     * position.
     * 
     * @param c
     * @return
     */
    @CompileStatic
    Matrix transform(Closure c) {
        
        if(c.maximumNumberOfParameters == 1) {
            return transformWithoutIndices(c)
        }
        else 
        if(c.maximumNumberOfParameters == 3) {
            return transformWithIndices(c)
        }
    }
    
    @CompileStatic
    Matrix transformWithoutIndices(Closure c) {
        final int rows = matrix.rowDimension
        final int cols = matrix.columnDimension
        double[][] newData = new double[rows][cols]
        for(int i=0; i<rows;++i) {
            double [] row = matrix.dataRef[i]
            for(int j=0; j<cols;++j) {
                newData[i][j] = (double)c(row[j])
            }
        }                    
        return new Matrix(new Array2DRowRealMatrix(newData))
    }
    
    @CompileStatic
    Matrix transformWithIndices(Closure c) {
        final int rows = matrix.rowDimension
        final int cols = matrix.columnDimension
        double[][] newData = new double[rows][cols]
        for(int i=0; i<rows;++i) {
            double [] row = matrix.dataRef[i]
            double [] newRow = newData[i]
            for(int j=0; j<cols;++j) {
                double value = row[j] // NOTE: embedding this direclty in call below causes VerifyError with CompileStatic
                newRow[j] = (double)c.call(value,i,j)
            }
        }                    
        return new Matrix(new Array2DRowRealMatrix(newData))
    }
    
    /**
     * Transform the given matrix by passing each row to the given
     * closure. If the closure accepts two arguments then the 
     * row index is passed as well. The closure must return a 
     * double array to replace the array passed in. If null
     * is returned then the data is left unchanged.
     * 
     * @param c Closure to process the data with.
     * @return  transformed Matrix
     */
    @CompileStatic
    Matrix transformRows(Closure c) {
        final int rows = matrix.rowDimension
        final int cols = matrix.columnDimension
        
        double[][] newData = new double[rows][cols]
        
        if(c.maximumNumberOfParameters == 1) {
            for(int i=0; i<rows;++i) {
                newData[i] = c(matrix.dataRef[i])
            }
        }
        else 
        if(c.maximumNumberOfParameters == 2) {
            for(int i=0; i<rows;++i) {
                newData[i] = c(matrix.dataRef[i], i)
            }
        }
        else
            throw new IllegalArgumentException("Closure must accept 1 or two arguments")
        
        return new Matrix(new Array2DRowRealMatrix(newData))
    }
    
    @CompileStatic
    void eachRow(Closure c) {
        if(c.maximumNumberOfParameters == 1) {
            for(double [] row in matrix.dataRef) {
                c(row)    
            }
        }
        else 
        if(c.maximumNumberOfParameters == 2) {
            int rowIndex = 0;
            for(double [] row in matrix.dataRef) {
                c(row, rowIndex)
                ++rowIndex
            }
        }
    }
    
    String toString() {
        matrix.data.collect { row -> 
            (row as List).join(",") 
        }.join("\n")
    }
}
