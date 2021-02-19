package gngs

import static org.junit.Assert.*

import org.junit.Test

class UtilsTest {

    List<Map> tableData = [
        [ person: 'Bob', age: 30 ],
        [ person: 'Tim', age: 50 ],
        [ person: 'Jane', age: 29 ],
        [ person: 'Fred', age: 45 ],
    ]
    
    StringWriter sw = new StringWriter()

    @Test
    public void 'basic table display'() {
        def result = Utils.table(tableData, out: sw).toString()
        println(result)
        assert result.startsWith('person')
        assert result.readLines()[0].endsWith('| age')
    }
    
    @Test
    public void 'table with top border'() {
        def result = Utils.table(tableData, topborder: true, out: sw).toString()
        println(result)
        assert result.startsWith('------')
        assert result.readLines()[1].endsWith('| age')
        assert result.readLines()[1].startsWith('person')
    }

    @Test
    public void 'table with top and side border'() {
        def result = Utils.table(tableData, topborder: true, border:true, out: sw).toString()
        println(result)
        assert result.startsWith(' ------')
        assert result.readLines()[1].endsWith('| age |')
        assert result.readLines()[1].startsWith('| person')
    }
 

    @Test
    void 'table with side border'() {
        
        def result = Utils.table(tableData, border: '|', out: sw).toString()
        
        println(result)

        assert result.startsWith('| person')
        assert result.readLines()[0].endsWith('| age |')
        assert result.readLines().contains('| Bob    | 30  |')
    }
}
