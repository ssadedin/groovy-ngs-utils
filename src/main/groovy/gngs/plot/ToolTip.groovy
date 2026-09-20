package gngs.plot

import java.awt.BasicStroke
import java.awt.Color
import java.awt.Font
import java.awt.FontMetrics
import java.awt.Graphics2D
import java.awt.RenderingHints
import java.awt.Stroke
import java.awt.geom.Ellipse2D
import java.awt.geom.Line2D
import java.awt.geom.Point2D
import java.awt.geom.Rectangle2D
import java.awt.geom.RoundRectangle2D

import de.erichseifert.gral.graphics.AbstractDrawable
import de.erichseifert.gral.graphics.Drawable
import de.erichseifert.gral.graphics.DrawingContext
import de.erichseifert.gral.plots.XYPlot
import de.erichseifert.gral.plots.axes.Axis
import de.erichseifert.gral.plots.axes.AxisRenderer
import de.erichseifert.gral.util.PointND

import groovy.transform.CompileStatic
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * Visual attributes controlling how tooltips are drawn.
 * <p>
 * A {@link Plot} carries a default instance which individual calls to
 * {@code showTooltips} may override on a per-series basis.
 *
 * @author Simon Sadedin
 */
@CompileStatic
class ToolTipStyle {

    /**
     * Explicit font. If null, the plot font is used, resized to {@link #fontSize}.
     */
    Font font = null

    float fontSize = 11.0f

    Color background = new Color(255, 255, 240, 240)

    Color border = new Color(110, 110, 110)

    Color textColor = Color.black

    Color leaderColor = new Color(80, 80, 80)

    /**
     * Horizontal padding between the text and the edge of the tooltip box
     */
    double padX = 7.0d

    /**
     * Vertical padding between the text and the edge of the tooltip box
     */
    double padY = 5.0d

    /**
     * Gap in pixels between the annotated data point and the nearest edge of
     * the tooltip box.
     */
    double distance = 26.0d

    double cornerRadius = 6.0d

    /**
     * Radius of the dot drawn over the annotated data point. Zero to omit it.
     */
    double markerRadius = 3.5d

    /**
     * Whether to draw a line connecting the data point to its tooltip. This is
     * what makes it viable for automatic placement to move a tooltip well away
     * from its point in order to avoid an overlap.
     */
    boolean leader = true

    /**
     * If set, text is wrapped at word boundaries so that the tooltip is no
     * wider than this, in pixels.
     */
    Double maxWidth = null

    /**
     * Create an independent copy of this style, with the given property values
     * applied over the top.
     *
     * @param overrides property values to set, which must all be properties of
     *                  {@code ToolTipStyle}
     */
    ToolTipStyle copy(Map overrides = null) {
        ToolTipStyle result = new ToolTipStyle(
            font: font,
            fontSize: fontSize,
            background: background,
            border: border,
            textColor: textColor,
            leaderColor: leaderColor,
            padX: padX,
            padY: padY,
            distance: distance,
            cornerRadius: cornerRadius,
            markerRadius: markerRadius,
            leader: leader,
            maxWidth: maxWidth
        )

        if(overrides) {
            // InvokerHelper.setProperties silently ignores properties that do
            // not exist, which would turn a typo into a mystery at render time
            for(Object key in overrides.keySet()) {
                if(result.hasProperty(key.toString()) == null)
                    throw new MissingPropertyException(
                        "Unknown tooltip style option '$key'. Valid options are: " +
                        styleOptions().join(', '), key.toString(), ToolTipStyle)
            }
            InvokerHelper.setProperties(result, overrides)
        }

        return result
    }

    /**
     * Names of the properties that may be given as style options
     */
    static List<String> styleOptions() {
        return new ToolTipStyle().properties.keySet()*.toString().grep { String name ->
            name != 'class'
        }.sort() as List<String>
    }
}

/**
 * Where a single tooltip should be placed relative to its data point.
 * <p>
 * A null {@link #angle} means the tooltip is placed automatically, by searching
 * for a position that does not collide with the data, the other tooltips, or
 * the edge of the plot.
 *
 * @author Simon Sadedin
 */
@CompileStatic
class ToolTipPlacement {

    /**
     * Direction from the data point to the tooltip, in degrees counter-clockwise
     * from east (so 90 is directly above the point). Null for automatic.
     */
    Double angle = null

    /**
     * Distance from the data point, overriding {@link ToolTipStyle#distance}
     */
    Double distance = null

    private static final Map<String, Double> ANCHOR_ANGLES = [
        'e'          : 0.0d,
        'east'       : 0.0d,
        'right'      : 0.0d,
        'ne'         : 45.0d,
        'north_east' : 45.0d,
        'n'          : 90.0d,
        'north'      : 90.0d,
        'above'      : 90.0d,
        'top'        : 90.0d,
        'nw'         : 135.0d,
        'north_west' : 135.0d,
        'w'          : 180.0d,
        'west'       : 180.0d,
        'left'       : 180.0d,
        'sw'         : 225.0d,
        'south_west' : 225.0d,
        's'          : 270.0d,
        'south'      : 270.0d,
        'below'      : 270.0d,
        'bottom'     : 270.0d,
        'se'         : 315.0d,
        'south_east' : 315.0d
    ]

    /**
     * Interpret a user supplied placement value.
     * <p>
     * Accepted forms are:
     * <ul>
     *   <li>null or {@code true} - place automatically</li>
     *   <li>a compass name such as {@code 'NW'}, or a direction such as
     *       {@code 'above'} / {@code 'left'}</li>
     *   <li>a number, interpreted as an angle in degrees counter-clockwise from east</li>
     *   <li>a map with {@code angle} and / or {@code distance} keys</li>
     * </ul>
     *
     * @throws IllegalArgumentException if the value cannot be interpreted
     */
    static ToolTipPlacement from(Object value) {

        if(value == null || value instanceof Boolean)
            return new ToolTipPlacement()

        if(value instanceof ToolTipPlacement)
            return (ToolTipPlacement)value

        if(value instanceof Number)
            return new ToolTipPlacement(angle: ((Number)value).doubleValue())

        if(value instanceof CharSequence) {
            String key = value.toString().toLowerCase().replaceAll(/[\s-]/, '_')
            Double angle = ANCHOR_ANGLES[key]
            if(angle == null)
                throw new IllegalArgumentException(
                    "Unknown tooltip anchor '$value'. Expected one of: " + ANCHOR_ANGLES.keySet().join(', '))
            return new ToolTipPlacement(angle: angle)
        }

        if(value instanceof Map) {
            Map map = (Map)value
            ToolTipPlacement result = new ToolTipPlacement()
            if(map.containsKey('anchor'))
                result.angle = from(map['anchor']).angle
            if(map['angle'] != null)
                result.angle = ((Number)map['angle']).doubleValue()
            if(map['distance'] != null)
                result.distance = ((Number)map['distance']).doubleValue()
            return result
        }

        throw new IllegalArgumentException(
            'Tooltip placement should be an anchor name, an angle, or a map of [angle: .., distance: ..], ' +
            'but was: ' + value.class.name)
    }
}

/**
 * A single tooltip that has been selected for display.
 *
 * @author Simon Sadedin
 */
@CompileStatic
class ToolTipAnnotation {

    /**
     * Index of the item within {@link ToolTipLayer#items} that this tooltip belongs to
     */
    int series

    /**
     * Index of the annotated point within its series
     */
    int index

    String text

    ToolTipPlacement placement = new ToolTipPlacement()

    ToolTipStyle style

    /**
     * View coordinates of the annotated data point, set when the tooltip is drawn
     */
    Point2D anchor = null

    /**
     * View bounds of the tooltip box, set when the tooltip is drawn. Useful for
     * understanding where automatic placement decided to put a tooltip.
     */
    Rectangle2D bounds = null
}

/**
 * Text that has been parsed, wrapped and measured, ready to draw.
 */
@CompileStatic
class LaidOutText {
    List<List<StyledRun>> lines
    double width
    double height
    double[] lineHeights
    int[] lineAscents
}

/**
 * Draws tooltips over an {@link XYPlot} as annotations on selected data points.
 * <p>
 * The layer is added to the plot as an ordinary GRAL component, but without a
 * layout constraint, which means GRAL will not attempt to position it and it is
 * drawn last, over the top of the plotted data.
 * <p>
 * Positions are resolved lazily, inside {@link #draw}, because the geometry of
 * a plot is not known until it has been laid out, which happens well after the
 * plot is constructed. Resolution follows exactly what GRAL itself does when it
 * renders a data point: the position of the plot area, offset by the position
 * the axis renderers assign to the data value.
 *
 * @author Simon Sadedin
 */
@CompileStatic
class ToolTipLayer extends AbstractDrawable {

    /**
     * Candidate directions tried when placing a tooltip automatically, in
     * preference order. The diagonals come first because they are least likely
     * to sit on top of a trend line.
     */
    private static final double[] AUTO_ANGLES = [
        45.0d, 135.0d, 315.0d, 225.0d,
        90.0d, 0.0d, 180.0d, 270.0d,
        22.5d, 67.5d, 112.5d, 157.5d, 202.5d, 247.5d, 292.5d, 337.5d
    ] as double[]

    /**
     * Multiples of the configured distance tried when placing automatically.
     * The leader line is what makes the larger values usable.
     */
    private static final double[] AUTO_DISTANCE_FACTORS = [1.0d, 1.9d, 2.8d] as double[]

    private static final double PENALTY_ESCAPE = 4.0d

    private static final double PENALTY_OVERLAP = 6.0d

    private static final double PENALTY_OBSTACLE = 900.0d

    XYPlot plot

    /**
     * Every XY item in the plot. Their data points and line segments are
     * treated as obstacles that automatic placement tries to avoid.
     */
    List<XYItem> items = []

    List<ToolTipAnnotation> annotations = []

    @Override
    void draw(DrawingContext context) {

        if(annotations.isEmpty())
            return

        Drawable plotArea = plot.getPlotArea()
        if(plotArea == null)
            return

        Rectangle2D plotBounds = plotArea.getBounds()
        if(plotBounds.width <= 0 || plotBounds.height <= 0)
            return

        // The axis renderers are what convert data values to view coordinates,
        // so without them there is nowhere to put a tooltip
        if(plot.getAxisRenderer(XYPlot.AXIS_X) == null || plot.getAxisRenderer(XYPlot.AXIS_Y) == null)
            return

        List<double[][]> positions = items.collect { XYItem item -> resolve(plotBounds, item) }

        // Line segments only exist for items that are actually drawn as lines
        List<double[]> segments = []
        items.eachWithIndex { XYItem item, int i ->
            if(!(item instanceof Lines) && !(item instanceof Line) && !(item instanceof Area))
                return

            double[][] pos = positions[i]
            double[] xs = pos[0]
            double[] ys = pos[1]
            for(int p = 1; p < xs.length; ++p) {
                segments.add([xs[p - 1], ys[p - 1], xs[p], ys[p]] as double[])
            }
        }

        List<Rectangle2D> fixed = fixedObstacles()
        List<Rectangle2D> placed = []

        Graphics2D g = context.getGraphics()
        Object antialiasing = g.getRenderingHint(RenderingHints.KEY_ANTIALIASING)
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON)
        Font fontOld = g.getFont()
        Color colorOld = g.getColor()
        Stroke strokeOld = g.getStroke()

        try {
            for(ToolTipAnnotation annotation in annotations) {

                double[][] pos = positions[annotation.series]
                if(annotation.index < 0 || annotation.index >= pos[0].length)
                    continue

                Point2D anchor = new Point2D.Double(pos[0][annotation.index], pos[1][annotation.index])

                ToolTipStyle style = annotation.style
                Font font = fontFor(style)
                LaidOutText text = layoutText(g, font, annotation.text, style)
                if(text.lines.isEmpty())
                    continue

                Rectangle2D box =
                    place(anchor, text, annotation, style, plotBounds, positions, segments, fixed, placed)

                placed.add(box)

                annotation.anchor = anchor
                annotation.bounds = box

                render(g, font, text, style, anchor, box)
            }
        }
        finally {
            g.setStroke(strokeOld)
            g.setColor(colorOld)
            g.setFont(fontOld)
            if(antialiasing != null)
                g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, antialiasing)
        }
    }

    /**
     * Convert the data values of an item to view coordinates, the same way that
     * {@code XYPlotArea2D} does when it renders them.
     *
     * @return a two element array holding the x and y view coordinates
     */
    private double[][] resolve(Rectangle2D plotBounds, XYItem item) {

        Axis axisX = plot.getAxis(XYPlot.AXIS_X)
        Axis axisY = plot.getAxis(XYPlot.AXIS_Y)
        AxisRenderer rendererX = plot.getAxisRenderer(XYPlot.AXIS_X)
        AxisRenderer rendererY = plot.getAxisRenderer(XYPlot.AXIS_Y)

        List xValues = item.x as List
        List yValues = item.y as List

        int n = Math.min(xValues.size(), yValues.size())
        double[] xs = new double[n]
        double[] ys = new double[n]

        for(int i = 0; i < n; ++i) {
            PointND<Double> px = rendererX.getPosition(axisX, (Number)xValues[i], true, false)
            PointND<Double> py = rendererY.getPosition(axisY, (Number)yValues[i], true, false)
            if(px == null || py == null) {
                // Off scale in a way the renderer cannot express: park it far
                // outside the plot so that it is never chosen or drawn
                xs[i] = Double.NaN
                ys[i] = Double.NaN
            }
            else {
                xs[i] = plotBounds.getMinX() + px.get(PointND.X)
                ys[i] = plotBounds.getMinY() + py.get(PointND.Y)
            }
        }

        return [xs, ys] as double[][]
    }

    /**
     * Parts of the plot other than the data that a tooltip should not cover
     */
    private List<Rectangle2D> fixedObstacles() {
        List<Rectangle2D> result = []

        if(plot.isLegendVisible() && plot.getLegend() != null)
            result.add(plot.getLegend().getBounds())

        if(plot.getTitle() != null)
            result.add(plot.getTitle().getBounds())

        return result.grep { Rectangle2D r -> r != null && r.width > 0 && r.height > 0 } as List<Rectangle2D>
    }

    private Font fontFor(ToolTipStyle style) {
        if(style.font != null)
            return style.font

        Font base = plot.getFont()
        return base == null ? new Font('SansSerif', Font.PLAIN, (int)style.fontSize) : base.deriveFont(style.fontSize)
    }

    private static Font fontFor(Font base, StyledRun run) {
        int style = Font.PLAIN
        if(run.bold)
            style |= Font.BOLD
        if(run.italic)
            style |= Font.ITALIC

        return style == Font.PLAIN ? base : base.deriveFont(style)
    }

    /**
     * Parse, optionally wrap, and measure the text of a tooltip.
     */
    private LaidOutText layoutText(Graphics2D g, Font baseFont, String source, ToolTipStyle style) {

        List<List<StyledRun>> lines = Markup.parse(source)

        if(style.maxWidth != null)
            lines = wrap(g, baseFont, lines, style.maxWidth - 2 * style.padX)

        int n = lines.size()
        double[] heights = new double[n]
        int[] ascents = new int[n]
        double width = 0
        double height = 0

        FontMetrics baseMetrics = g.getFontMetrics(baseFont)

        for(int i = 0; i < n; ++i) {
            List<StyledRun> line = lines[i]
            double lineWidth = 0
            double lineHeight = 0
            int ascent = 0
            for(StyledRun run in line) {
                FontMetrics fm = g.getFontMetrics(fontFor(baseFont, run))
                lineWidth += fm.stringWidth(run.text)
                lineHeight = Math.max(lineHeight, (double)fm.getHeight())
                ascent = Math.max(ascent, fm.getAscent())
            }

            if(line.isEmpty()) {
                lineHeight = (double)baseMetrics.getHeight()
                ascent = baseMetrics.getAscent()
            }

            heights[i] = lineHeight
            ascents[i] = ascent
            width = Math.max(width, lineWidth)
            height += lineHeight
        }

        return new LaidOutText(
            lines: lines,
            width: width + 2 * style.padX,
            height: height + 2 * style.padY,
            lineHeights: heights,
            lineAscents: ascents
        )
    }

    /**
     * Re-flow lines so that no line exceeds the given width, breaking at spaces.
     */
    private List<List<StyledRun>> wrap(Graphics2D g, Font baseFont, List<List<StyledRun>> lines, double maxWidth) {

        List<List<StyledRun>> result = []

        for(List<StyledRun> line in lines) {

            List<StyledRun> current = []
            double currentWidth = 0

            for(StyledRun run in line) {
                FontMetrics fm = g.getFontMetrics(fontFor(baseFont, run))

                // Words carry their trailing whitespace so that it is preserved
                // when the word stays on the line, and discarded when it does not
                for(String word in run.text.split(/(?<=\s)(?=\S)/)) {

                    double wordWidth = fm.stringWidth(word)

                    if(currentWidth > 0 && (currentWidth + wordWidth) > maxWidth) {
                        result.add(current)
                        current = []
                        currentWidth = 0
                        word = word.replaceAll(/^\s+/, '')
                        wordWidth = fm.stringWidth(word)
                    }

                    StyledRun last = current.isEmpty() ? null : current[-1]
                    if(last != null && last.bold == run.bold && last.italic == run.italic)
                        last.text = last.text + word
                    else
                        current.add(new StyledRun(text: word, bold: run.bold, italic: run.italic))

                    currentWidth += wordWidth
                }
            }

            result.add(current)
        }

        return result
    }

    /**
     * Decide where the tooltip box goes.
     * <p>
     * An explicitly specified angle is honoured as given. Otherwise candidate
     * positions are scored and the least bad one chosen.
     */
    private Rectangle2D place(Point2D anchor, LaidOutText text, ToolTipAnnotation annotation,
                              ToolTipStyle style, Rectangle2D plotBounds, List<double[][]> positions,
                              List<double[]> segments, List<Rectangle2D> fixed, List<Rectangle2D> placed) {

        double distance = annotation.placement.distance != null ? annotation.placement.distance : style.distance

        if(annotation.placement.angle != null)
            return boxAt(anchor, text, annotation.placement.angle, distance)

        // Only obstacles that could possibly be reached by this tooltip matter.
        // Filtering them once per tooltip keeps placement viable on plots with
        // very large numbers of points.
        double reach = distance * AUTO_DISTANCE_FACTORS[AUTO_DISTANCE_FACTORS.length - 1] +
                       Math.max(text.width, text.height)
        Rectangle2D searchArea = new Rectangle2D.Double(
            anchor.x - reach, anchor.y - reach, 2 * reach, 2 * reach)

        List<double[]> nearPoints = []
        for(double[][] pos in positions) {
            double[] xs = pos[0]
            double[] ys = pos[1]
            for(int i = 0; i < xs.length; ++i) {
                if(!Double.isNaN(xs[i]) && searchArea.contains(xs[i], ys[i]))
                    nearPoints.add([xs[i], ys[i]] as double[])
            }
        }

        List<double[]> nearSegments = segments.grep { double[] s ->
            !Double.isNaN(s[0]) && !Double.isNaN(s[2]) && searchArea.intersectsLine(s[0], s[1], s[2], s[3])
        } as List<double[]>

        Rectangle2D best = null
        double bestPenalty = Double.MAX_VALUE

        for(double factor in AUTO_DISTANCE_FACTORS) {
            for(double angle in AUTO_ANGLES) {
                Rectangle2D candidate = boxAt(anchor, text, angle, distance * factor)

                double penalty = penalty(candidate, plotBounds, nearPoints, nearSegments, fixed, placed, style) +
                                 (factor - 1.0d) * distance * 0.4d

                if(penalty < bestPenalty) {
                    bestPenalty = penalty
                    best = candidate
                }

                if(bestPenalty == 0.0d)
                    return best
            }
        }

        return best
    }

    /**
     * The box of the given size, offset from the anchor in the given direction
     * such that its nearest edge is approximately {@code distance} away.
     */
    private static Rectangle2D boxAt(Point2D anchor, LaidOutText text, double angleDegrees, double distance) {
        double radians = Math.toRadians(angleDegrees)
        double centreX = anchor.x + Math.cos(radians) * (distance + text.width / 2.0d)

        // Screen y runs downwards, but angles are expressed the way a reader
        // expects, ie: 90 degrees is above the point
        double centreY = anchor.y - Math.sin(radians) * (distance + text.height / 2.0d)

        return new Rectangle2D.Double(
            centreX - text.width / 2.0d, centreY - text.height / 2.0d, text.width, text.height)
    }

    /**
     * Score a candidate position. Lower is better, zero is unobstructed.
     */
    private double penalty(Rectangle2D candidate, Rectangle2D plotBounds, List<double[]> points,
                           List<double[]> segments, List<Rectangle2D> fixed, List<Rectangle2D> placed,
                           ToolTipStyle style) {

        double result = 0

        Rectangle2D inside = candidate.createIntersection(plotBounds)
        double escaped = area(candidate) - area(inside)
        result += escaped * PENALTY_ESCAPE

        for(Rectangle2D other in placed) {
            result += area(other.createIntersection(candidate)) * PENALTY_OVERLAP
        }

        for(Rectangle2D other in fixed) {
            result += area(other.createIntersection(candidate)) * PENALTY_OVERLAP
        }

        double margin = Math.max(2.0d, style.markerRadius)
        Rectangle2D expanded = new Rectangle2D.Double(
            candidate.x - margin, candidate.y - margin,
            candidate.width + 2 * margin, candidate.height + 2 * margin)

        for(double[] point in points) {
            if(expanded.contains(point[0], point[1]))
                result += PENALTY_OBSTACLE
        }

        for(double[] segment in segments) {
            if(candidate.intersectsLine(segment[0], segment[1], segment[2], segment[3]))
                result += PENALTY_OBSTACLE
        }

        return result
    }

    private static double area(Rectangle2D r) {
        return (r.width <= 0 || r.height <= 0) ? 0.0d : r.width * r.height
    }

    private void render(Graphics2D g, Font baseFont, LaidOutText text, ToolTipStyle style,
                        Point2D anchor, Rectangle2D box) {

        if(style.leader) {
            g.setColor(style.leaderColor)
            g.setStroke(new BasicStroke(1.0f))
            g.draw(new Line2D.Double(anchor, attachPoint(anchor, box)))
        }

        if(style.markerRadius > 0) {
            g.setColor(style.leaderColor)
            double r = style.markerRadius
            g.fill(new Ellipse2D.Double(anchor.x - r, anchor.y - r, 2 * r, 2 * r))
        }

        RoundRectangle2D outline = new RoundRectangle2D.Double(
            box.x, box.y, box.width, box.height, style.cornerRadius, style.cornerRadius)

        if(style.background != null) {
            g.setColor(style.background)
            g.fill(outline)
        }

        if(style.border != null) {
            g.setColor(style.border)
            g.setStroke(new BasicStroke(1.0f))
            g.draw(outline)
        }

        g.setColor(style.textColor)

        double y = box.y + style.padY
        text.lines.eachWithIndex { List<StyledRun> line, int i ->
            double x = box.x + style.padX
            for(StyledRun run in line) {
                Font font = fontFor(baseFont, run)
                g.setFont(font)
                g.drawString(run.text, (float)x, (float)(y + text.lineAscents[i]))
                x += g.getFontMetrics(font).stringWidth(run.text)
            }
            y += text.lineHeights[i]
        }
    }

    /**
     * The point on the box outline nearest the data point, so that the leader
     * line stops at the edge of the box rather than running underneath it.
     */
    private static Point2D attachPoint(Point2D anchor, Rectangle2D box) {
        double x = Math.max(box.minX, Math.min(anchor.x, box.maxX))
        double y = Math.max(box.minY, Math.min(anchor.y, box.maxY))

        boolean inside = x > box.minX && x < box.maxX && y > box.minY && y < box.maxY
        if(inside)
            return new Point2D.Double(box.centerX, box.centerY)

        return new Point2D.Double(x, y)
    }
}
