package gngs.plot

import org.junit.Test

class MarkupTest {

    @Test
    void 'plain text is a single line with a single unstyled run'() {
        List<List<StyledRun>> lines = Markup.parse('hello world')

        assert lines.size() == 1
        assert lines[0].size() == 1
        assert lines[0][0].text == 'hello world'
        assert !lines[0][0].bold
        assert !lines[0][0].italic
    }

    @Test
    void 'null and empty input produce no lines'() {
        assert Markup.parse(null) == []
        assert Markup.parse('') == []
    }

    @Test
    void 'bold and italic tags create styled runs'() {
        List<List<StyledRun>> lines = Markup.parse('a<b>b</b><i>c</i>d')

        assert lines.size() == 1
        assert lines[0]*.text == ['a', 'b', 'c', 'd']
        assert lines[0]*.bold == [false, true, false, false]
        assert lines[0]*.italic == [false, false, true, false]
    }

    @Test
    void 'strong and em are synonyms for bold and italic'() {
        assert Markup.parse('<strong>x</strong>')[0][0].bold
        assert Markup.parse('<em>x</em>')[0][0].italic
    }

    @Test
    void 'styles nest'() {
        List<StyledRun> runs = Markup.parse('<b>a<i>b</i>c</b>')[0]

        assert runs*.text == ['a', 'b', 'c']
        assert runs*.bold == [true, true, true]
        assert runs*.italic == [false, true, false]
    }

    @Test
    void 'br breaks lines in all its spellings'() {
        assert Markup.parse('a<br>b').size() == 2
        assert Markup.parse('a<br/>b').size() == 2
        assert Markup.parse('a<br />b').size() == 2
        assert Markup.parse('a<BR>b').size() == 2
    }

    @Test
    void 'a literal newline breaks the line'() {
        List<List<StyledRun>> lines = Markup.parse('a\nb')

        assert lines.size() == 2
        assert lines[0][0].text == 'a'
        assert lines[1][0].text == 'b'
    }

    @Test
    void 'paragraphs become separate lines without leading blanks'() {
        List<List<StyledRun>> lines = Markup.parse('<p>one</p><p>two</p>')

        assert lines.size() == 2
        assert lines[0][0].text == 'one'
        assert lines[1][0].text == 'two'
    }

    @Test
    void 'entities are decoded'() {
        assert Markup.parse('a &amp; b')[0]*.text.join('') == 'a & b'
        assert Markup.parse('&lt;tag&gt;')[0]*.text.join('') == '<tag>'
        assert Markup.parse('&quot;q&quot;')[0]*.text.join('') == '"q"'
        assert Markup.parse('a&nbsp;b')[0]*.text.join('') == 'a b'
    }

    @Test
    void 'an unknown entity is left alone'() {
        assert Markup.parse('100 &copyright; x')[0]*.text.join('') == '100 &copyright; x'
        assert Markup.parse('a & b')[0]*.text.join('') == 'a & b'
    }

    @Test
    void 'unknown tags are stripped but their content is kept'() {
        assert Markup.parse('<span style="color:red">kept</span>')[0][0].text == 'kept'
        assert Markup.parse('<div><h1>title</h1></div>')[0][0].text == 'title'
    }

    @Test
    void 'a stray closing tag is ignored rather than failing'() {
        List<StyledRun> runs = Markup.parse('a</b>b')[0]

        assert runs*.text.join('') == 'ab'
        assert runs.every { !it.bold }
    }

    @Test
    void 'an unterminated tag is treated as literal text'() {
        assert Markup.parse('a <b')[0]*.text.join('') == 'a <b'
    }

    @Test
    void 'whitespace inside text is preserved'() {
        // Tooltips often align values with spaces, so unlike HTML we must not collapse them
        assert Markup.parse('a    b')[0][0].text == 'a    b'
    }

    @Test
    void 'strip returns plain text with newlines for breaks'() {
        assert Markup.strip('<b>NA12878</b><br>depth: 62x') == 'NA12878\ndepth: 62x'
        assert Markup.strip('plain') == 'plain'
        assert Markup.strip(null) == ''
    }

    @Test
    void 'a realistic tooltip parses as expected'() {
        List<List<StyledRun>> lines = Markup.parse('<b>NA12878</b><br>depth: <i>62x</i><br>chr1:1,000,000')

        assert lines.size() == 3
        assert lines[0][0].text == 'NA12878'
        assert lines[0][0].bold
        assert lines[1]*.text == ['depth: ', '62x']
        assert lines[1][1].italic
        assert lines[2][0].text == 'chr1:1,000,000'
    }
}
