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

object bayesFit extends App {

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
  val data : Seq[(Double, Double)] = boneLengths.zip(statures)
  
    

  // We need this line to seed breeze's random number generator
  implicit val randBasisBreeze: breeze.stats.distributions.RandBasis = rng.breezeRandBasis
  
  // prior specification
  val dSlope = Gaussian(3, 1)
  val dIntercept = Gaussian(1700, 100)
  val dSigma = Uniform(1, 100)
  println(boneLengthsMean)

  //val a = 0.2
  //val b = 3
  //val sigma2 = 0.5
  //val errorDist = breeze.stats.distributions.Gaussian(0, sigma2)(rng.breezeRandBasis)
  //val data = for (x <- 0 until 100) yield {
    //(x.toDouble, a * x + b + errorDist.draw())
  //}

  case class Parameters(a: Double, b: Double, sigma2: Double)

  implicit object tuple3ParameterConversion
      extends ParameterConversion[Tuple3[Double, Double, Double], Parameters] {
    def from(p: Parameters): Tuple3[Double, Double, Double] = (p.a, p.b, p.sigma2)
    def to(t: Tuple3[Double, Double, Double]): Parameters =
      Parameters(a = t._1, b = t._2, sigma2 = t._3)
  }
  case class LikelihoodEvaluator(data: Seq[(Double, Double)])
      extends MHDistributionEvaluator[Parameters] {

    override def logValue(theta: MHSample[Parameters]): Double = {

      val likelihoods = for ((x, y) <- data) yield {
        val likelihood = breeze.stats.distributions.Gaussian(
          theta.parameters.a * (x - boneLengthsMean) + theta.parameters.b,
          theta.parameters.sigma2
        )

        likelihood.logPdf(y)
      }
      likelihoods.sum
    }
  }

  object PriorEvaluator extends MHDistributionEvaluator[Parameters] {

    val priorDistA = dSlope
    val priorDistB = dIntercept
    val priorDistSigma = dSigma

    override def logValue(theta: MHSample[Parameters]): Double = {
      priorDistA.logPdf(theta.parameters.a)
      +priorDistB.logPdf(theta.parameters.b)
      +priorDistSigma.logPdf(theta.parameters.sigma2)
    }
  }
  val posteriorEvaluator = ProductEvaluator(PriorEvaluator, LikelihoodEvaluator(data))

  val genA = GaussianRandomWalkProposal(0.05, "rw-a-0.1").forType[Double]
  val genB = GaussianRandomWalkProposal(30, "rw-b-0.5").forType[Double]
  val genSigma = GaussianRandomWalkProposal(0.1, "rw-sigma-0.01").forType[Double]

  val parameterGenerator = MHProductProposal(genA, genB, genSigma).forType[Parameters]

  val identProposal = MHIdentityProposal.forType[Double]
  val noiseOnlyGenerator =
    MHProductProposal(identProposal, identProposal, genSigma).forType[Parameters]

  val mixtureGenerator = MHMixtureProposal((0.1, noiseOnlyGenerator), (0.9, parameterGenerator))
  val chain = MetropolisHastings(mixtureGenerator, posteriorEvaluator)

  val logger = MHSampleLogger[Parameters]()
  val initialSample = MHSample(Parameters(3.0, 1700.0, 10.0), generatedBy = "initial")

  val mhIterator = chain.iterator(initialSample, logger)
  val samples = mhIterator.drop(1000).take(5000).toIndexedSeq
  val meanAndVarianceA = meanAndVariance(samples.map(_.parameters.a))
  println(
    s"Estimates for parameter a: mean = ${meanAndVarianceA.mean}, var = ${meanAndVarianceA.variance}"
  )
  val meanAndVarianceB = meanAndVariance(samples.map(_.parameters.b))
  println(
    s"Estimates for parameter b: mean = ${meanAndVarianceB.mean}, var = ${meanAndVarianceB.variance}"
  )
  val meanAndVarianceSigma2 = meanAndVariance(samples.map(_.parameters.sigma2))
  println(
    s"Estimates for parameter sigma2: mean = ${meanAndVarianceSigma2.mean}, var = ${meanAndVarianceSigma2.variance}"
  )
  println("Acceptance ratios: " + logger.samples.acceptanceRatios)
  println(
    "acceptance ratios over the last 100 samples: " + logger.samples.takeLast(100).acceptanceRatios
  )
}
