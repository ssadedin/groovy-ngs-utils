/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2014 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
import groovy.transform.CompileStatic;

import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.RealMatrix;

/**
 * A proxy object representing a column in a matrix. 
 * <p>
 * The data in a {@link Matrix} is stored natively in row format. That is,
 * each row is stored as a native Java array of double values. This makes
 * accessing data by row very efficient, but doesn't give you an easy way to
 * pass around or treat a column of values as a collection without 
 * first copying them to another data structure. This class wraps 
 * an {@link Iterable} interface around a column of values without actually copying
 * the data. It does this keeps a reference to the underlying matrix and 
 * implements iteration and random access (via square bracket notation) 
 * by reflecting values into the appropriate column of the underlying
 * Matrix.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
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
        sourceMatrix.matrix.rowDimension
    }
    
    Object asType(Class c) {
        if(c == List) {
            return sourceMatrix.matrix.getColumn(columnIndex) as List
        }
        else
        if(c == double[]) {
            return sourceMatrix.matrix.getColumn(columnIndex)
        }
        else {
            return super.asType(c)
        }
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
    
    String toString() {
        "[" + this.collect {it}.join(",") + "]"
    }
}

/**
 * Wraps an Apache-Commons-Math matrix of double values with a 
 * Groovy interface for convenient access. Because it wraps the
 * underlying matrix as a delegate, all the original methods of
 * the Commons-Math implementation are available directly, along with
 * Groovy-enhanced versions.
 * <p>
 * The most basic enhancements come in the form of random access operators 
 * that allowthe Matrix class to be referenced using square-bracket notation:
 * <pre>
 * Matrix m = new Matrix(2,2,[1,2,3,4])
 * assert m[0][0] == 2
 * assert m[1][1] == 4
 * </pre>
 * The rows of the Matrix are directly accessible simply by using
 * square-bracket indexing:
 * <pre>
 * assert m[0] == [1,2]
 * assert m[1] == [3,4]
 * </pre>
 * The columns are accessed by using an empty first index:
 * <pre>
 * assert m[][0] == [1,3]
 * assert m[][1] == [2,4]
 * </pre>
 * Rows and columns can both be treated as normal Groovy collections:
 * <pre>
 * assert m[0].collect { it.any { it > 2 }  } == [ false, true ]
 * assert m[][0].collect { it > 1 } == [ false, true ]
 * </pre>
 * Note that in the above code, both row-wise and column-wise access
 * occurs without copying any data.
 * <p>
 * Transforming the whole matrix can be done using <code>transform</code>:
 * <pre>
 * assert m.transform { it * 2 } == Matrix(2,2,[2,4,6,8])
 * </pre>
 * As an option, row and column indexes are available as well:
 * <pre>
 * assert m.transform { value, row, column -> value * 2 } == Matrix(2,2,[2,4,6,8])
 * </pre>
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Matrix {
    
    static { 
        double[][].metaClass.toMatrix = { new Matrix(delegate) }
        
        def originalMethod = double[][].metaClass.getMetaMethod("asType", Class)
        double[][].metaClass.asType = { arg -> arg == Matrix.class ? delegate.toMatrix() : originalMethod(arg)}
    }
    
    /**
     * How many rows are displayed in toString() and other calls that format output
     */
    static final int DISPLAY_ROWS = 50
    
    @Delegate
    Array2DRowRealMatrix matrix

    public Matrix(int rows, int columns) {
        matrix = new Array2DRowRealMatrix(rows, columns)
    }
    
    public Matrix(MatrixColumn... columns) {
        matrix = new Array2DRowRealMatrix(columns[0].size(), columns.size())
        
    }
    
    public Matrix(double [][] values) {
        matrix = new Array2DRowRealMatrix(values, false)
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
        if(n instanceof Number)
            return matrix.dataRef[(int)n]
        else
        if(n instanceof List) {
            List<Number> l = (List)n
            if(l.size() == 0) // Seems to happen with m[][2] type syntax
                return getColumns()
            else {
                return subsetRows(l)
            }
        }
        else
        if(n instanceof Iterable) {
            return subsetRows((n))
        }
        else
        if(n.class.isArray()) {
            return subsetRows(n as Collection<Number>)
        }
        else {
            throw new IllegalArgumentException("Cannot subset rows by type: " + n?.class?.name)
        }
    }
    
    @CompileStatic
    Object subsetRows(Iterable<Number> i) {
        List indices = new ArrayList(this.matrix.rowDimension)
        i.each { Number n -> indices.add(n.toInteger()) }
        return matrix.dataRef[indices]
    }
    
    @CompileStatic
    void putAt(int n, Object values) {
       matrix.dataRef[n] = (values as double[])
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
    
    
    /**
     * Filter the rows of this matrix and return 
     * a Matrix as a result
     * 
     * @param   c   a Closure to evaluate
     * 
     * @return  Matrix for which the closure c returns a non-false value
     */
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
     * @param c A closure taking 1 or 3 arguments (data value, or data value, row,
     *          column)
     * @return  A matrix reflecting the transformed data values
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
    private Matrix transformWithoutIndices(Closure c) {
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
    private Matrix transformWithIndices(Closure c) {
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
    
    /**
     * Shorthand to give a familiar function to R users
     */
    List<Long> which(Closure c) {
        this.findIndexValues(c)
    }
    
    String toString() {
        if(matrix.rowDimension<DISPLAY_ROWS) {
            int rowCount = 0
            return "${matrix.rowDimension}x${matrix.columnDimension} Matrix:\n"+ matrix.data.collect { row -> 
                (rowCount++) + ":\t" + (row as List).join(",\t") 
            }.join("\n")
        }
        else {
            int omitted = matrix.rowDimension-DISPLAY_ROWS
            int rowCount = 0
            return "${matrix.rowDimension}x${matrix.columnDimension} Matrix:\n"+ matrix.data[0..DISPLAY_ROWS/2].collect { row -> 
                ((rowCount++) + ":").padRight(6) + (row as List).join(",\t") 
            }.join("\n") + "\n... ${omitted} rows omitted ...\n" + matrix.data[-(DISPLAY_ROWS/2)..-1].collect { row -> 
                (((rowCount++)+omitted-1) + ":").padRight(6) + (row as List).join(",\t") 
            }.join("\n")
        }
    }
}
