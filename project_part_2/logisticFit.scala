//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using dep "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.sampling.MHSample
import scalismo.sampling.MHDistributionEvaluator
import scalismo.sampling.evaluators.ProductEvaluator
import scalismo.sampling.MHProposalGenerator
import scalismo.sampling.proposals.GaussianRandomWalkProposal
import scalismo.sampling.proposals.MHProductProposal
import scalismo.sampling.ParameterConversion
import scalismo.sampling.algorithms.MetropolisHastings
import scalismo.sampling.loggers.MHSampleLogger
import scalismo.sampling.proposals.MHMixtureProposal
import scalismo.sampling.proposals.MHIdentityProposal
import breeze.stats.meanAndVariance
import scalismo.plot.data.DataFrame
import breeze.stats.distributions.Gaussian
import breeze.stats.distributions.LogNormal
import breeze.stats.distributions.Uniform
import breeze.macros.expand.sequence
import breeze.stats.distributions.Bernoulli
import scala.math._

object logisticFit extends App {

  scalismo.initialize()
  implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

  // Read data using scalismo-plot
  val dataFrame = DataFrame.fromCSV(java.io.File("Project data/statureAndSex.csv")).get

  val boneLengths : Seq[Double] = dataFrame.column("bone-length")
    .values
    .map(cellValue => cellValue.asContinuous.value)

  val statures : Seq[Double] = dataFrame.column("stature")
    .values
    .map(cellValue => cellValue.asContinuous.value)

  val trochanterDistances : Seq[Double] = dataFrame.column("trochanter-distance")
    .values
    .map(cellValue => cellValue.asContinuous.value)

  val sexs : Seq[String] = dataFrame.column("sex")
    .values
    .map(cellValue => cellValue.toString)

  val ids : Seq[String] = dataFrame.column("id")
    .values
    .map(cellValue => cellValue.toString)


  val boneLengthsMean = boneLengths.sum / boneLengths.length
  val staturesMean = statures.sum / statures.length
  val demeanedBoneLengths = boneLengths.map(_ - boneLengthsMean)
  val demeanedStatures = statures.map(_ - staturesMean)
  //val data : Seq[(Double, Double)] = demeanedBoneLengths.zip(demeanedStatures)
  val data : Seq[(Double, String)] = boneLengths.zip(sexs)
  
  
  // the inverse logit function
  def invLogit(x: Double): Double = 1.0 / (1.0 + math.exp(-x))

  // We need this line to seed breeze's random number generator
  implicit val randBasisBreeze: breeze.stats.distributions.RandBasis = rng.breezeRandBasis
  
  // prior specification
  val dSlope = Gaussian(0.5, 0.1)
  val dIntercept = Gaussian(0, 0.1)


  //val a = 0.2
  //val b = 3
  //val sigma2 = 0.5
  //val errorDist = breeze.stats.distributions.Gaussian(0, sigma2)(rng.breezeRandBasis)
  //val data = for (x <- 0 until 100) yield {
    //(x.toDouble, a * x + b + errorDist.draw())
  //}

  case class Parameters(a: Double, b: Double)

  implicit object tuple3ParameterConversion
      extends ParameterConversion[Tuple2[Double, Double], Parameters] {
    def from(p: Parameters): Tuple2[Double, Double] = (p.a, p.b)
    def to(t: Tuple2[Double, Double]): Parameters =
      Parameters(a = t._1, b = t._2)
  }
  case class LikelihoodEvaluator(data: Seq[(Double, String)])
      extends MHDistributionEvaluator[Parameters] {

    override def logValue(theta: MHSample[Parameters]): Double = {

      val likelihoods = for ((x, y) <- data) yield {
        var i = 0
        if (y.equals("m")) {
          
          val likelihood = invLogit(theta.parameters.a * (x - boneLengthsMean) + theta.parameters.b)
          log(likelihood)
           
        } else {
          
          val likelihood = 1 - invLogit(theta.parameters.a * (x - boneLengthsMean) + theta.parameters.b)
          log(likelihood)

        }
        

      }
      likelihoods.sum
    }
  }

  object PriorEvaluator extends MHDistributionEvaluator[Parameters] {

    val priorDistA = dSlope
    val priorDistB = dIntercept

    override def logValue(theta: MHSample[Parameters]): Double = {
      priorDistA.logPdf(theta.parameters.a)
      +priorDistB.logPdf(theta.parameters.b)
    }
  }
  val posteriorEvaluator = ProductEvaluator(PriorEvaluator, LikelihoodEvaluator(data))

  val genA = GaussianRandomWalkProposal(0.1, "rw-a-0.1").forType[Double]
  val genB = GaussianRandomWalkProposal(0.1, "rw-b-0.5").forType[Double]

  val parameterGenerator = MHProductProposal(genA, genB).forType[Parameters]

  val identProposal = MHIdentityProposal.forType[Double]
  val noiseOnlyGenerator =
    MHProductProposal(identProposal, identProposal).forType[Parameters]

  val mixtureGenerator = MHMixtureProposal((0.01, noiseOnlyGenerator), (0.99, parameterGenerator))
  val chain = MetropolisHastings(mixtureGenerator, posteriorEvaluator)

  val logger = MHSampleLogger[Parameters]()
  val initialSample = MHSample(Parameters(0.5, 0), generatedBy = "initial")

  val mhIterator = chain.iterator(initialSample, logger)
  val samples = mhIterator.drop(1000).take(10000).toIndexedSeq
  val meanAndVarianceA = meanAndVariance(samples.map(_.parameters.a))
  println(
    s"Estimates for parameter a: mean = ${meanAndVarianceA.mean}, var = ${meanAndVarianceA.variance}"
  )
  val meanAndVarianceB = meanAndVariance(samples.map(_.parameters.b))
  println(
    s"Estimates for parameter b: mean = ${meanAndVarianceB.mean}, var = ${meanAndVarianceB.variance}"
  )
  println("Acceptance ratios: " + logger.samples.acceptanceRatios)
  println(
    "acceptance ratios over the last 100 samples: " + logger.samples.takeLast(100).acceptanceRatios
  )
}
