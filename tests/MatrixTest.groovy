import static org.junit.Assert.*;

import org.junit.Test;


class MatrixTest {

    @Test
    public void test() {
        
        def m = new Matrix(3,4)
        
        m[2][3] = 9
        
        assert m[2][3] == 9
        
        assert m[2,3] == 9
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
}
