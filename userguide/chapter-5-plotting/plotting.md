Visualisation of data is an extremely common task in working with genomic data. Unfortunately,
the libraries available for generating plots in Java are mostly geared towards
applications and are not as easy to use or feature rich for use with
scientific data.

To ease plot generation, Groovy NGS includes the [Gral](https://github.com/eseifert/gral)
plotting library to support generation of plots by default, allowing any
APIs from this library to be called to  create visualisations. As with
many Java based plotting APIs however, it is not well integrated with data types
typically used in data analysis, requiring many API calls to generate a simple
plot. Groovy NGS therefore includes a set of wrappers that make generation
of useful plots directly from raw data very ergonomic. These APIs are based
around the plotting libraries supported in 
the [BeakerX JupyterLab Extenions](https://github.com/twosigma/beakerx), allowing
the same plots that are generated in notebooks to then be easily used offline in command 
line applications.

# Basic Form of Plots

Plots are created by creating a `Plot` object:

```groovy
import gngs.plot.*

p = new Plot()
```

Most attributes of the plots are passed using Groovy's [Map Constructor](https://docs.groovy-lang.org/latest/html/documentation/#_named_parameters)
idiom, such that plot-level features are mostly configured in the constructor by passing
them as attributes. An example of a plot with title, axis labels and X and Y ranges set for the 
plot area is as follows:

```groovy
p = new Plot(
    title: 'An example plot',
    xLabel: 'This is the x axis',
    yLabel: 'This is the y axis',
    xBound: [0, 10],
    yBound: [0, 30]
)
```

# Saving Plots

Once you have created a plot, you will likely want to save it to a file.

This can easily be achieved by calling the `save` function:

```groovy
p.save('example_plot.png')
```

By default, images are saved resolution of 1024 x 800, but you can set this
explicitly by passing `width` and `height` properties:

```groovy
p.save('example_plot.png', width: 1600, height: 1200)
```

At this stage, these will save an empty plot, because we have not added any data to it.


# Adding Data to Plots

Data is added to plots by creating data objects and using Groovy's `<<` operator on the plot 
object. The data objects themselves are created using the same map-constructor style syntax,
and can be `Line`, `Bars`, `Points`, or `Area` objects:

```groovy
p << new Line(
    x: [1,2,3,4,5],
    y: [1,4,9,16,25]
)
```

:include-image: example_plot.png {scale: 0.5}

Note that in most of these data objects, the actual position of the data points is expressed
as `x` and `y` attributes of the objects.

# Specifying Colors

Further attributes can be specified such as the color of the objects. The color can be
specified as [java.awt.Color](https://docs.oracle.com/javase%2F7%2Fdocs%2Fapi%2F%2F/java/awt/Color.html) 
objects, for which the built in predefined colors are often useful:

```groovy
import java.awt.Color

p << new Line(
    x: [1,2,3,4,5],
    y: [1,4,9,16,25],
    color: Color.red
)

p << new Line(
    x: [1,2,3,4,5],
    y: [1,8,27,64,125],
    color: Color.blue
)
```

# Adding Legends

To add a legend to the plot, specify the `displayName` attribute on each
data object to be drawn:

```groovy
import java.awt.Color

p << new Line(
    x: [1,2,3,4,5],
    y: [1,4,9,16,25],
    color: Color.red,
    displayName: 'Square of Numbers'
)

p << new Line(
    x: [1,2,3,4,5],
    y: [1,8,27,64,125],
    color: Color.blue,
    displayName: 'Cube of Numbers'
)
```

:include-image: example_plot2.png {scale: 0.5}

# Distribution Plots

As they are very common to show overall data properties, Groovy NGS supports
two kinds of distribution plots that are especially easy to create: Density plots
and Histograms.

```groovy
p = new Plot(title: 'Guassian Distribution')

p << new gngs.plot.bx.Density.Area(
       data: (1..100).collect { 5 +  r.nextGaussian() }
    )
```

:include-image: density.png {scale: 0.5}
    
    
```groovy
h = new Histogram(
   title: "Gaussian Distribution",
   data: (1..500).collect { 5 +  r.nextGaussian() }
)
```

:include-image: histogram.png {scale: 0.5}
    
# Accessing Generated Images     

If you want to add more to the generated plots or manipulate the image in
additional ways, you can generate the [BufferedImage](https://docs.oracle.com/javase/8/docs/api/java/awt/image/BufferedImage.html) 
directly and use standard Java features to draw on it:

```groovy
BufferedImage image = p.image
Graphics2D imageGraphics = image.createGraphics()
imageGraphics.drawRect(10,10, 20,20)
...
```

As with saving the plot, you can specify the resolution:

```groovy
BufferedImage image = p.getImage(1600, 1200)
...
```






