import groovy.transform.CompileStatic;

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.RealMatrix;

class MatrixColumn {
    
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
                return index<max
            }
            
            Double next() {
                return sourceMatrix.dataRef[index++][columnIndex]
            }
            
            void remove() { throw new UnsupportedOperationException() }
        }
    }
}

class Matrix {
    
    @Delegate
    Array2DRowRealMatrix matrix

    public Matrix(int rows, int columns) {
        matrix = new Array2DRowRealMatrix(rows, columns)
    }
    
    public Matrix(Array2DRowRealMatrix m) {
        matrix = m
    }
     
    @CompileStatic
    MatrixColumn col(int n) {
//        matrix.getColumn(n)
        new MatrixColumn(columnIndex:n, sourceMatrix: this)
    }
    
    @CompileStatic
    double [] row(int n) {
        matrix.getRow(n)
    }
    
    @CompileStatic
    Object getAt(Object n) {
        if(n instanceof Integer)
            return matrix.dataRef[(int)n]
        else
        if(n instanceof List) {
            List l = (List)n
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
