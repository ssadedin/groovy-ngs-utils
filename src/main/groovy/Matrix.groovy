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
	
	MatrixColumn() {
		name = "C"+columnIndex
	}
	
	String name
    
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
    
    boolean equals(Object o) {
        int i = 0
        return (o.every { it == this[i++] })
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
		
//		println "Setting Matrix meta class properties ...."
        double[][].metaClass.toMatrix = { new Matrix(delegate) }
        
        def originalMethod = double[][].metaClass.getMetaMethod("asType", Class)
        double[][].metaClass.asType = { arg -> arg == Matrix.class ? delegate.toMatrix() : originalMethod(arg)}
        
        def original2 = Array2DRowRealMatrix.metaClass.getMetaMethod("asType", Class)
        Array2DRowRealMatrix.metaClass.asType = { arg -> arg == Matrix.class ? new Matrix(arg) : original2(arg)}
        
        def originalMultiply = Integer.metaClass.getMetaMethod("multiply", Class)
        Integer.metaClass.multiply = { arg -> arg instanceof Matrix ? arg.multiply(delegate) : originalMultiply(arg)}
    }
	
    /**
     * How many rows are displayed in toString() and other calls that format output
     */
    static final int DISPLAY_ROWS = 50
    
    @Delegate
    Array2DRowRealMatrix matrix
	
	List<String> names = []
	
    public Matrix(int rows, int columns) {
        matrix = new Array2DRowRealMatrix(rows, columns)
    }
    
    @CompileStatic
    public Matrix(MatrixColumn... columns) {
        int rows = columns[0].size()
        double[][] newData =  new double[rows][]
        for(int i=0; i<rows;++i) {
            newData[i] = (double[])columns.collect { MatrixColumn c -> c[i] } as double[]
        }
        matrix = new Array2DRowRealMatrix(newData)
        this.names = columns.collect { MatrixColumn c -> c.name }
    }
    
    public Matrix(List<Iterable> rows) {
        double [][] data = null
        int rowCount = 0
        for(r in rows) {
            double[] rowData = r.collect{ (double)it }
            if(data == null)
                data = new double[rows.size()][rowData.size()]
            data[rowCount++] = rowData
        }
        matrix = new Array2DRowRealMatrix(data, false)
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
        new MatrixColumn(columnIndex:n, sourceMatrix: this, name: names[n])
    }
    
    List<MatrixColumn> getColumns(List<String> names) {
        names.collect { this.names.indexOf(it) }.collect { int index ->
             assert index >= 0; col(index) 
        }
    }
    
    List<MatrixColumn> getColumns() {
        (0..<matrix.columnDimension).collect { col(it) }
    }
    
    @CompileStatic
    double [] row(int n) {
        matrix.getRow(n)
    }
    
	/**
	 * Implementation of the [] operator. Adds several different behaviors:
	 * <ul>
	 *    <li>Plain old indexing returns a row: <code>m[4]</code> returns 5th row of matrix.
	 *    <li>Double indexing returns a cell: <code>m[4][5]</code> returns 6th column of 4th row.
	 *    <li>Empty index returns a column: <code>m[][6]</code> returns 7th column
	 *    <li>List (or any iterable) index returns rows matching indices:
	 *    </ul>
	 * <pre>
	 * Matrix m = new Matrix([1..80], 10, 8)
	 * m[2..4] == [ [ 9..16 ], [17..24], [25..32] ]
	 * @param n
	 * @return
	 */
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
    
	/**
	 * Return a subset of the rows indicated by the indices in the given iterable
	 * (Note that the indices don't need to be consecutive or monotonic).
	 * @param i
	 * @return
	 */
    @CompileStatic
    Object subsetRows(Iterable<Number> i) {
        List<Integer> indices = new ArrayList(this.matrix.rowDimension)
        i.each { Number n -> indices.add(n.toInteger()) }
		
		double [][] result = new double[indices.size()][this.matrix.columnDimension]
		int destRowIndex = 0
		for(int srcRowIndex in indices) {
			System.arraycopy(this.matrix.dataRef[srcRowIndex], 0, result[destRowIndex++], 0, this.matrix.columnDimension)
		}
        return result
    }
    
    @CompileStatic
    void putAt(int n, Object values) {
       matrix.dataRef[n] = (values as double[])
    }
    
	/**
	 * Specialization of <code>collect</code>: if 1 arg then
	 * just pass the row, if 2 args, pass the row AND the index.
	 * 
	 * @param c	Closure to execute for each row in the matrix
	 * @return	results collected
	 */
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
    
    Matrix multiply(double d) {
        new Matrix(this.matrix.scalarMultiply(d))
    }
    
    Matrix multiply(Matrix m) {
        new Matrix(this.matrix.preMultiply(m.dataRef))
    }
    
    Matrix plus(Matrix m) {
        new Matrix(this.matrix.add(m.matrix))
    }
      
    Matrix minus(Matrix m) {
        new Matrix(this.matrix.subtract(m.matrix))
    }
    
     Matrix divide(double d) {
        new Matrix(this.matrix.scalarMultiply(1/d))
    }
    
    Matrix plus(double x) {
        new Matrix(this.matrix.scalarAdd(x))
    }
    
    Matrix minus(double x) {
        new Matrix(this.matrix.scalarMinus(x))
    }
    
    Matrix transpose() {
        new Matrix(this.matrix.transpose())
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
