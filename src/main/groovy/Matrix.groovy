import graxxia.MatrixColumn;

import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

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

 /**
  * This class aliases the Graxxia Matrix class into the default namespace
  * for backwards compatibility.
  * 
  * @author simon
  */
class Matrix extends graxxia.Matrix {

    public Matrix(Array2DRowRealMatrix m) {
        super(m);
    }

    public Matrix(double[]... values) {
        super(values);
    }

    public Matrix(int rows, int columns, double... matrixData) {
        super(rows, columns, matrixData);
    }

    public Matrix(int rows, int columns, List<Double> data) {
        super(rows, columns, data);
    }

    public Matrix(int rows, int columns) {
        super(rows, columns);
    }

    public Matrix(Iterable<Iterable> rows) {
        super(rows);
    }

    public Matrix(MatrixColumn... sourceColumns) {
        super(sourceColumns);
    }
}