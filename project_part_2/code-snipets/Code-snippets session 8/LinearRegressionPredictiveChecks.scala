//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using dep "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.plot.data.DataFrame
import scalismo.plot.data.DataRow
import scalismo.plot.data.DataFrame.Column
import scalismo.plot.plottarget.PlotTargets.plotTargetBrowser
import breeze.stats.distributions.Gaussian
import breeze.stats.distributions.LogNormal
import breeze.stats.distributions.Rand.FixedSeed.randBasis


object LinearRegressionPredictiveChecks {

  case class Sample(length : Double, stature : Double, slope : Double, intercept : Double, sigma : Double)

  def main(args: Array[String]): Unit = {

    // prior specification
    val dSlope = Gaussian(1, 1)
    val dIntercept = Gaussian(0, 2)
    val dSigma = LogNormal(0, 0.25)

    // samples from the prior
    val samples = for length <- 300 until 500 yield 
      val slope = dSlope.sample()
      val intercept = dIntercept.sample()
      val sigma = dSigma.sample()
      val stature = slope * length + intercept + Gaussian(0, sigma).sample()
      Sample(length, stature, slope, intercept, sigma)

    // Create a scatterplot and a histogram from the samples
    val df = DataFrame.fromColumns(Seq(
        Column.ofContinuous(samples.map(_.length), "length"),
        Column.ofContinuous(samples.map(_.stature), "stature")
    ))
    df.plot.scatterPlot("length", "stature", title="stature from length").show()
    df.plot.histogram("stature", title="stature").show()


    // create a line plot from the first sample
    val firstSample = samples.head
    val lengths = Seq.range(300, 500).map(_.toDouble)
    val stature = lengths.map(length => firstSample.slope * length + firstSample.intercept)
    DataFrame.fromColumns(Seq(
        Column.ofContinuous(lengths, "length"),
        Column.ofContinuous(stature, "stature")
    )).plot.linePlot("length", "stature", title="stature from length").show()  
  }

    
}
