package gngs.plot

import groovy.transform.CompileStatic

/**
 * A run of text within a tooltip that shares a single font style.
 * <p>
 * Runs are the output of {@link Markup#parse(String)}: a tooltip is parsed into
 * a list of lines, and each line is a list of runs.
 *
 * @author Simon Sadedin
 */
@CompileStatic
class StyledRun {

    String text

    boolean bold

    boolean italic

    String toString() {
        String style = (bold ? 'b' : '') + (italic ? 'i' : '')
        return style ? "$text($style)" : text
    }
}

/**
 * Parser for the very small subset of HTML that is commonly used in BeakerX
 * tooltips.
 * <p>
 * This is deliberately <i>not</i> a HTML parser. It recognises a fixed set of
 * inline tags and entities, and silently discards everything else it does not
 * understand, so that markup carried over from a notebook degrades to clean
 * text rather than appearing literally in a rendered image.
 * <p>
 * Supported:
 * <ul>
 *   <li>{@code <b>}, {@code <strong>} - bold</li>
 *   <li>{@code <i>}, {@code <em>} - italic</li>
 *   <li>{@code <br>}, {@code <br/>} - line break</li>
 *   <li>{@code <p>}, {@code </p>} - line break (a paragraph boundary)</li>
 *   <li>the entities {@code &lt; &gt; &amp; &quot; &apos; &nbsp;}</li>
 *   <li>literal newline characters, which break lines just like {@code <br>}</li>
 * </ul>
 * Any other tag is stripped, leaving its text content in place. Stray closing
 * tags are ignored rather than treated as an error, on the basis that a
 * malformed tooltip should not prevent a plot from rendering.
 * <p>
 * Note that, unlike HTML, whitespace inside text is preserved exactly. Tooltips
 * frequently contain deliberately aligned values, and collapsing runs of spaces
 * would destroy that.
 * <p>
 * Example:
 * <pre>
 * Markup.parse('&lt;b&gt;NA12878&lt;/b&gt;&lt;br&gt;depth: 62x')
 * // =&gt; [[NA12878(b)], [depth: 62x]]
 * </pre>
 *
 * @author Simon Sadedin
 */
@CompileStatic
class Markup {

    private static final Map<String, String> ENTITIES = [
        'lt'   : '<',
        'gt'   : '>',
        'amp'  : '&',
        'quot' : '"',
        'apos' : "'",
        'nbsp' : ' '
    ]

    private static final char LT = '<' as char

    private static final char AMP = '&' as char

    private static final char NEWLINE = '\n' as char

    /**
     * Parse the given source text into lines of styled runs.
     *
     * @param src   tooltip text, which may contain simple markup, or null
     * @return      a list of lines, each being a list of {@link StyledRun}. Empty
     *              if the source is null or contains no renderable text.
     */
    static List<List<StyledRun>> parse(String src) {

        List<List<StyledRun>> lines = []

        if(src == null)
            return lines

        List<StyledRun> line = []
        StringBuilder buf = new StringBuilder()
        int bold = 0
        int italic = 0

        int i = 0
        final int n = src.length()
        while(i < n) {
            char c = src.charAt(i)

            if(c == LT) {
                int close = src.indexOf('>', i)
                if(close < 0) {
                    // Unterminated tag: there is nothing sensible to do except
                    // treat the remainder as literal text
                    buf.append(src.substring(i))
                    break
                }

                String tag = src.substring(i + 1, close).trim().toLowerCase()
                boolean closing = tag.startsWith('/')
                String name = tagName(closing ? tag.substring(1) : tag)

                if(name == 'b' || name == 'strong') {
                    flush(line, buf, bold, italic)
                    bold = closing ? Math.max(0, bold - 1) : bold + 1
                }
                else
                if(name == 'i' || name == 'em') {
                    flush(line, buf, bold, italic)
                    italic = closing ? Math.max(0, italic - 1) : italic + 1
                }
                else
                if(name == 'br') {
                    flush(line, buf, bold, italic)
                    lines.add(line)
                    line = []
                }
                else
                if(name == 'p') {
                    // A paragraph boundary only breaks the line if there is
                    // something on the current line to break away from
                    if(!line.isEmpty() || buf.length() > 0) {
                        flush(line, buf, bold, italic)
                        lines.add(line)
                        line = []
                    }
                }
                // else: unknown tag, strip it

                i = close + 1
            }
            else
            if(c == AMP) {
                int semi = src.indexOf(';', i)
                String entity = null
                if(semi > i && (semi - i) <= 8)
                    entity = src.substring(i + 1, semi).toLowerCase()

                String replacement = entity == null ? null : ENTITIES[entity]
                if(replacement != null) {
                    buf.append(replacement)
                    i = semi + 1
                }
                else {
                    buf.append(c)
                    ++i
                }
            }
            else
            if(c == NEWLINE) {
                flush(line, buf, bold, italic)
                lines.add(line)
                line = []
                ++i
            }
            else {
                buf.append(c)
                ++i
            }
        }

        flush(line, buf, bold, italic)
        if(!line.isEmpty())
            lines.add(line)

        return lines
    }

    /**
     * Return the given text with all markup removed, as a single plain string
     * with line breaks rendered as newline characters.
     */
    static String strip(String src) {
        return parse(src).collect { List<StyledRun> line ->
            line*.text.join('')
        }.join('\n')
    }

    /**
     * Extract the tag name from the body of a tag, ie: everything up to the
     * first whitespace or self closing slash.
     */
    private static String tagName(String tag) {
        int end = tag.length()
        for(int i = 0; i < tag.length(); ++i) {
            char c = tag.charAt(i)
            if(Character.isWhitespace(c) || c == ('/' as char)) {
                end = i
                break
            }
        }
        return tag.substring(0, end)
    }

    /**
     * Move any buffered text into the current line as a run with the given style.
     */
    private static void flush(List<StyledRun> line, StringBuilder buf, int bold, int italic) {
        if(buf.length() == 0)
            return

        line.add(new StyledRun(text: buf.toString(), bold: bold > 0, italic: italic > 0))
        buf.setLength(0)
    }
}
