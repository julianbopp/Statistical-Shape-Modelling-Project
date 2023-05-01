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

object LogisticRegression {

  // the inverse logit function
  def invLogit(x: Double): Double = 1.0 / (1.0 + math.exp(-x))

  enum Sex: 
    case Male extends Sex
    case Female extends Sex
  
  
  case class Sample(length : Double, p : Double, sex : Sex, slope : Double, intercept : Double)


  def main(args: Array[String]): Unit = {
    import breeze.stats.distributions.Rand.FixedSeed.randBasis

    // prior specification
    val dSlope = Gaussian(3, 1)
    val dIntercept = Gaussian(-1000, 100)

    // samples from the prior
    val samples = for length <- 300 until 500 yield 
      val slope = dSlope.sample()
      val intercept = dIntercept.sample()
      val p = invLogit(slope * (length) + intercept)
      val sex = breeze.stats.distributions.Bernoulli(p).sample()
      Sample(length, p, if sex == true then Sex.Male else Sex.Female, slope, intercept)


    val df = DataFrame.fromColumns(Seq(
        Column.ofContinuous(samples.map(sample => sample.length), "length"),
        Column.ofContinuous(samples.map(sample => sample.p), "p"),
        Column.ofNominals(samples.map(sample => sample.sex.toString()), "sex")
    ))
    df.plot.scatterPlot("length", "p", title="logit", colorField = "sex").show()
  }

    
}
