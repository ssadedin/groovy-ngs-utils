import static org.junit.Assert.*;

import org.junit.Test;


class MatrixTest {

    @Test
    public void test() {
        
        def m = new Matrix(3,4)
        
        m[2][3] = 9
        
        assert m[2][3] == 9
    }
    
    @Test
    void testColumn() {
        def m = new Matrix(4,2, [0d,1d,
                                 2d,3d,
                                 4d,5d,
                                 6d,6d])
        
        assert m.col(0)[0] == 0d
        assert m.col(0)[1] == 2d
        
        println m.col(0).grep { println(it); it > 0 }
        
        println m.col(0).find { it > 0 }
        
        println m.col(0).findIndexValues { it > 0 }
        
        println m.col(0)[0,2]
        
        
        println m.col(1).collect { it %2 }
    }
    
    @Test
    void testColumnMean() {
        def m = new Matrix(4,2, [0d,1d,
                                 2d,3d,
                                 4d,5d,
                                 6d,6d])
       
        assert Stats.mean(m.col(0)) == 3
        assert Stats.mean(m.col(1)) == 15 / 4
        
        assert Stats.mean(m[][0]) == 3
     }
    
    @Test
    void testColumnAccess() {
        def m = new Matrix(4,3, [0d,1d,3d,
                                 2d,3d,4d,
                                 4d,5d,6d,
                                 6d,6d,6d])
        
        assert m.columns.size() == 3
        
        println "m[][2] = " + m.columns[2]
        println "m[][2] = " + m[][2]
        println "m[][2][3] = " + m[][2][3]
        assert m[][2][3] == 6d
    }
    
    @Test
    void bigMatrix() {
        Matrix m = new Matrix(60,2)
        for(int i=0; i<m.rowDimension;++i) {
            m[i] = [i,i+1]
        }
        
        println "Matrix is $m"
        
        assert m.toString().contains("10 rows omitted")
    }
    
    @Test
    void testListAccess() {
        def m = new Matrix(4,3, [0d,1d,3d,
                                 2d,3d,4d,
                                 4d,5d,6d,
                                 6d,6d,6d])
        def sub = m[ [1,2]]         
        assert sub[0][0] == 2
        assert sub[1][2] == 6
        
        int [] indices = [1,2] as int[]
        sub = m[indices]
        assert sub[0][0] == 2
        assert sub[1][2] == 6
        
        long [] lndices = [1,2] as long[]
        sub = m[lndices]
        assert sub[0][0] == 2
        assert sub[1][2] == 6
         
        List lndices2 = [1L,2L]
        sub = m[lndices2]
        assert sub[0][0] == 2
        assert sub[1][2] == 6
    }
    
    @Test
    void testTypeConversion() {
        double [][] x = [
         [2d,3d,4d] as double[],
         [4d,5d,6d] as double[]
        ]
        
        Matrix m = x.toMatrix()
        assert m[0] == [2d,3d,4d]
        
        Matrix m2 = x as Matrix
        assert m2[0] == [2d,3d,4d]
        
        List l = x as List
        assert l[0] == [2d,3d,4d]
        
        def c = m[][1]
        
        assert c as List == [3d,5d]
    }
    
}
