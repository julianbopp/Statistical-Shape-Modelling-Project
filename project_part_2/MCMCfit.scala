//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using dep "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"
//
import scalismo.io.StatisticalModelIO
import scalismo.io.LandmarkIO
import scalismo.ui.api.ScalismoUI
import scalismo.geometry._
import scalismo.common.PointId
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.common.UnstructuredPointsDomain
import scalismo.common.interpolation.NearestNeighborInterpolator3D
import scalismo.common.UnstructuredPointsDomain1D
import scalismo.common.UnstructuredPointsDomain3D
import scalismo.statisticalmodel.PointDistributionModel
import scalismo.statisticalmodel.MultivariateNormalDistribution

import scalismo.transformations._

import scalismo.sampling._
import scalismo.sampling.proposals._
import scalismo.sampling.parameters._
import scalismo.sampling.evaluators._
import scalismo.sampling.loggers.MHSampleLogger
import scalismo.sampling.algorithms.MetropolisHastings

import breeze.linalg.DenseVector
import breeze.linalg.DenseMatrix

import scalismo.io.{MeshIO, StatisticalModelIO, LandmarkIO}
import scalismo.common._
import scalismo.mesh._
import scalismo.numerics.UniformMeshSampler3D
import scalismo.io.{MeshIO, StatisticalModelIO, LandmarkIO}

import scalismo.ui.api._
import breeze.stats.distributions.MultivariateGaussian
import breeze.linalg._
import scalismo.plot.data.*
import scalismo.plot.plottarget.PlotTargets.plotTargetBrowser

object MCMCfit extends App {

  implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)
  implicit val randBasisBreeze: breeze.stats.distributions.RandBasis = rng.breezeRandBasis
  scalismo.initialize()

  val ui = ScalismoUI()

  // Load and display shape model
  val model =
    StatisticalModelIO.readStatisticalTriangleMeshModel3D(new java.io.File("ssm.h5")).get
  val modelGroup = ui.createGroup("model")
  val modelView = ui.show(modelGroup, model, "model")
  modelView.referenceView.opacity = 0.5

  // Load and display fragments
  val targetGroup = ui.createGroup("fragments")
  val targetMeshes : Seq[TriangleMesh[_3D]] = (0 until 10).flatMap { i =>
    val filename = s"fragment-$i.stl" 
    val mesh = MeshIO.readMesh(new java.io.File(s"Project data/fragments/$filename")).get
    //ui.show(targetGroup, mesh, s"fragment-$i")
    Some(mesh)
  }

  def computeCenterOfMass(mesh: TriangleMesh[_3D]): Point[_3D] = {
    val normFactor = 1.0 / mesh.pointSet.numberOfPoints
    mesh.pointSet.points.foldLeft(Point(0, 0, 0))((sum, point) => sum + point.toVector * normFactor)
  }

  val rotationCenter = computeCenterOfMass(model.mean)

  case class Parameters(poseAndShapeParameters: PoseAndShapeParameters, noiseStddev: Double)
  object Parameters {

    implicit object parameterConversion
        extends ParameterConversion[
          Tuple2[PoseAndShapeParameters, Double],
          Parameters
        ] {
      def from(p: Parameters): Tuple2[PoseAndShapeParameters, Double] =
        (p.poseAndShapeParameters, p.noiseStddev)
      def to(t: Tuple2[PoseAndShapeParameters, Double]): Parameters =
        Parameters(t._1, t._2)
    }

    def poseTransformationForParameters(
        translationParameters: TranslationParameters,
        rotationParameters: RotationParameters,
        centerOfRotation: Point[_3D]
    ): TranslationAfterRotation[_3D] = {
      TranslationAfterRotation3D(
        Translation3D(translationParameters.translationVector),
        Rotation3D(rotationParameters.angles, centerOfRotation)
      )
    }
  }

  case class PriorEvaluator(model: PointDistributionModel[_3D, TriangleMesh])
      extends MHDistributionEvaluator[Parameters] {

    val translationPrior = breeze.stats.distributions.Gaussian(0.0, 5.0)
    val rotationPrior = breeze.stats.distributions.Gaussian(0, 0.1)
    val noisePrior = breeze.stats.distributions.LogNormal(1, 0.25)

    override def logValue(sample: MHSample[Parameters]): Double = {
      val poseAndShapeParameters = sample.parameters.poseAndShapeParameters
      val translationParameters = poseAndShapeParameters.translationParameters
      val rotationParameters = poseAndShapeParameters.rotationParameters

      model.gp.logpdf(poseAndShapeParameters.shapeParameters.coefficients) +
        translationPrior.logPdf(translationParameters.translationVector.x) +
        translationPrior.logPdf(translationParameters.translationVector.y) +
        translationPrior.logPdf(translationParameters.translationVector.z) +
        rotationPrior.logPdf(rotationParameters.angles._1) +
        rotationPrior.logPdf(rotationParameters.angles._2) +
        rotationPrior.logPdf(rotationParameters.angles._3) +
        noisePrior.logPdf(sample.parameters.noiseStddev)
    }
  }

  case class CorrespondenceEvaluator(
      model: PointDistributionModel[_3D, TriangleMesh],
      rotationCenter: Point[_3D],
      targetMesh: TriangleMesh[_3D]
  ) extends MHDistributionEvaluator[Parameters] {

    // we extract the points and build a model from only the points
    //val (refPoints, targetPoints) = correspondences.unzip

    //val newDomain = UnstructuredPointsDomain3D(refPoints.toIndexedSeq)
    //val modelOnLandmarkPoints = model.newReference(newDomain, NearestNeighborInterpolator3D())

    val mesh : TriangleMesh[_3D] = model.reference
    val downsampledMesh = mesh.operations.decimate(targetedNumberOfVertices = 250)
    val downsampledShapeModel = model.newReference(
      newReference = downsampledMesh,
      interpolator = TriangleMeshInterpolator3D()
    )

    val downsampledTargetMesh = targetMesh.operations.decimate(targetedNumberOfVertices = 250)
    override def logValue(sample: MHSample[Parameters]): Double = {

      val poseTransformation = Parameters.poseTransformationForParameters(
        sample.parameters.poseAndShapeParameters.translationParameters,
        sample.parameters.poseAndShapeParameters.rotationParameters,
        rotationCenter
      )
      val modelCoefficients = sample.parameters.poseAndShapeParameters.shapeParameters.coefficients
      val noiseModel = MultivariateGaussian(DenseVector.zeros[Double](3), DenseMatrix.eye[Double](3) * Math.abs(sample.parameters.noiseStddev))
      val transformedModel = downsampledShapeModel
        .instance(sample.parameters.poseAndShapeParameters.shapeParameters.coefficients)
        .transform(poseTransformation)

      val likelihoods = for (pointOnTarget <- downsampledTargetMesh.pointSet.points) yield {
        val u = transformedModel.operations.closestPointOnSurface(pointOnTarget).point - pointOnTarget
        noiseModel.logPdf(u.toBreezeVector)
      }

      val loglikelihood = likelihoods.sum
      loglikelihood
    }
  }

  val likelihoodEvaluator = CorrespondenceEvaluator(model, rotationCenter, targetMeshes(9))
  val priorEvaluator = PriorEvaluator(model).cached

  val posteriorEvaluator = ProductEvaluator(priorEvaluator, likelihoodEvaluator)

  val rotationProposal = MHProductProposal(
    GaussianRandomWalkProposal(0.001, "rx").forType[Double],
    GaussianRandomWalkProposal(0.001, "ry").forType[Double],
    GaussianRandomWalkProposal(0.01, "rz").forType[Double]
  ).forType[RotationParameters]

  val translationProposal = MHProductProposal(
    GaussianRandomWalkProposal(0.2, "tx").forType[Double],
    GaussianRandomWalkProposal(0.2, "ty").forType[Double],
    GaussianRandomWalkProposal(0.2, "tz").forType[Double]
  ).forType[TranslationParameters]

  val shapeProposalLeading =
    GaussianRandomWalkProposal(0.07, "shape-0-5")
      .partial(0 until 5)
      .forType[ShapeParameters]
  val shapeProposalRemaining =
    GaussianRandomWalkProposal(0.2, "shape-6-")
      .partial(6 until model.rank)
      .forType[ShapeParameters]

  val identTranslationProposal =
    MHIdentityProposal.forType[TranslationParameters]
  val identRotationProposal = MHIdentityProposal.forType[RotationParameters]
  val identShapeProposal = MHIdentityProposal.forType[ShapeParameters]

  val poseAndShapeTranslationOnlyProposal =
    MHProductProposal(
      translationProposal,
      identRotationProposal,
      identShapeProposal
    )
      .forType[PoseAndShapeParameters]
      .relabel("translation-only")
  val poseAndShapeRotationOnlyProposal =
    MHProductProposal(
      identTranslationProposal,
      rotationProposal,
      identShapeProposal
    )
      .forType[PoseAndShapeParameters]
      .relabel("rotation-only")
  val poseAndShapeLeadingShapeOnlyProposal =
    MHProductProposal(
      identTranslationProposal,
      identRotationProposal,
      shapeProposalLeading
    )
      .forType[PoseAndShapeParameters]
      .relabel("shape-leading-only")

  val poseAndShapeRemainingShapeOnlyProposal =
    MHProductProposal(
      identTranslationProposal,
      identRotationProposal,
      shapeProposalRemaining
    )
      .forType[PoseAndShapeParameters]
      .relabel("shape-trailing-only")

  val mixturePoseAndShapeProposal = MHMixtureProposal(
    (0.2, poseAndShapeTranslationOnlyProposal),
    (0.2, poseAndShapeRotationOnlyProposal),
    (0.3, poseAndShapeLeadingShapeOnlyProposal),
    (0.3, poseAndShapeRemainingShapeOnlyProposal)
  )
  val noiseProposal = GaussianRandomWalkProposal(0.1, "noise").forType[Double]
  val identNoiseProposal = MHIdentityProposal.forType[Double]
  val identPoseAndShapeProposal =
    MHIdentityProposal.forType[PoseAndShapeParameters]

  val noiseOnlyProposal =
    MHProductProposal(identPoseAndShapeProposal, noiseProposal)
      .forType[Parameters]
  val poseAndShapeOnlyProposal =
    MHProductProposal(mixturePoseAndShapeProposal, identNoiseProposal)
      .forType[Parameters]
  val fullproposal = MHMixtureProposal(
    (0.9, poseAndShapeOnlyProposal),
    (0.1, noiseOnlyProposal)
  )

  val logger = MHSampleLogger[Parameters]()
  val chain = MetropolisHastings(fullproposal, posteriorEvaluator)

  val initialParameters = Parameters(
    PoseAndShapeParameters(
      TranslationParameters(EuclideanVector3D(0, 0, 0)),
      RotationParameters((0.0, 0.0, 0.0)),
      ShapeParameters(DenseVector.zeros[Double](model.rank))
    ),
    noiseStddev = 10
  )

  val mhIterator = chain.iterator(MHSample(initialParameters, "inital"), logger)

  val samplingIterator = for ((sample, iteration) <- mhIterator.zipWithIndex) yield {
    println("iteration " + iteration)
    if (iteration % 500 == 0) {
      val poseAndShapeParameters = sample.parameters.poseAndShapeParameters
      val poseTransformation = Parameters.poseTransformationForParameters(
        poseAndShapeParameters.translationParameters,
        poseAndShapeParameters.rotationParameters,
        rotationCenter
      )
      modelView.shapeModelTransformationView.shapeTransformationView.coefficients =
        poseAndShapeParameters.shapeParameters.coefficients
      modelView.shapeModelTransformationView.poseTransformationView.transformation =
        poseTransformation

    }
    sample
  }

  val numOfIters = 100 
  val samples = samplingIterator.drop(0).take(numOfIters).toIndexedSeq
  println(logger.samples.acceptanceRatios)

  val logValues = samples.map(sample => posteriorEvaluator.logValue(sample))
  //println(logValues)

  val bestSample = samples.maxBy(posteriorEvaluator.logValue)

  val bestPoseAndShapeParameters = bestSample.parameters.poseAndShapeParameters
  val bestPoseTransformation = Parameters.poseTransformationForParameters(
    bestPoseAndShapeParameters.translationParameters,
    bestPoseAndShapeParameters.rotationParameters,
    rotationCenter
  )

  val bestFit = model
    .instance(bestPoseAndShapeParameters.shapeParameters.coefficients)
    .transform(bestPoseTransformation)
  val resultGroup = ui.createGroup("result")

  ui.show(resultGroup, bestFit, "best fit")

  val iterationsInt = 1 to numOfIters
  val iterations = iterationsInt.map(x => x.toDouble)
  val df = DataFrame(
    Seq(
      DataFrame.Column.ofContinuous(iterations, "iterations"),
      DataFrame.Column.ofContinuous(logValues, "logValues"),
    )
  )

  df.plot.linePlot(
    x = "iterations",
    y = "logValues",
    title = "test"
  ).show()

  
    val fileIds = 0 until 10

    def invLogit(x: Double): Double = 1.0 / (1.0 + math.exp(-x))
    // Computing the measurements
    val measurements = for (fileId <- fileIds) yield {
      println(s"processing $fileId")

      val landmarkFileA =
        new java.io.File( s"Project data/fragments/fits/fit${fileId}a.json")
      val landmarkFileB =
        new java.io.File( s"Project data/fragments/fits/fit${fileId}b.json")

      val landmarkA = LandmarkIO.readLandmarksJson3D(landmarkFileA).get
      val landmarkB = LandmarkIO.readLandmarksJson3D(landmarkFileB).get
      
      val LA = landmarkA.find(lm => lm.id == s"fit${fileId}a").get
      val LB = landmarkB.find(lm => lm.id == s"fit${fileId}b").get


      val length = (LA.point - LB.point).norm
      println("sex")
      val p = invLogit(0.495 * (length- 431) + 0.0338)
      println(p)
    }
  
}
