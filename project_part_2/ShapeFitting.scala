import scalismo.io.MeshIO
import scalismo.io.StatisticalModelIO
import scalismo.ui.api.*
import scalismo.sampling.ParameterConversion
import scalismo.sampling.parameters.*
import scalismo.sampling.*
import scalismo.sampling.evaluators.*
import scalismo.sampling.proposals.*
import scalismo.geometry.*
import scalismo.transformations.*
import scalismo.statisticalmodel.PointDistributionModel
import scalismo.mesh.TriangleMesh
import scalismo.sampling.algorithms.MetropolisHastings
import scalismo.common.interpolation.*
import breeze.linalg.DenseVector
import breeze.linalg.DenseMatrix
import breeze.stats.distributions.Rand.FixedSeed.randBasis
import scalismo.sampling.loggers.MHSampleLogger
import scalismo.plot.data.DataFrame
import scalismo.plot.data.DataFrame.Column
import scalismo.plot.plottarget.PlotTargets.plotTargetBrowser

object ShapeFitting {

    
  implicit val rng : scalismo.utils.Random = scalismo.utils.Random(42)
  val ui = ScalismoUI()


  //--------------------------------------------------------------------------------------------
  // Definition of parameters 
  //--------------------------------------------------------------------------------------------

  // The PoseAndShapeParameters are part of Scalismo
  case class Parameters(
      poseAndShapeParameters: PoseAndShapeParameters,
      noiseStddev: Double
  )

  // One more conversion needed as noise is additional parameter not available in Scalismo
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

    // Convenience function for getting the transformations from
    // the parameters
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

  //--------------------------------------------------------------------------------------------
  // Priors 
  //--------------------------------------------------------------------------------------------

  // We lump the priors all in one evaluator. We could, of course also have several 
  // separate evaluators, which we combine using a ProductEvaluator
  case class PriorEvaluator(model: PointDistributionModel[_3D, TriangleMesh])
      extends MHDistributionEvaluator[Parameters] {
      
    val translationPrior = breeze.stats.distributions.Gaussian(0.0, 5000.0) // applies for x,y and z
    val rotationPrior = breeze.stats.distributions.Uniform(-Math.PI, Math.PI)// applies to all rotation angles
    val noisePrior = breeze.stats.distributions.LogNormal(0, 0.5)

    override def logValue(sample: MHSample[Parameters]): Double = {
      val poseAndShapeParameters = sample.parameters.poseAndShapeParameters
      val translationParameters = poseAndShapeParameters.translationParameters
      val rotationParameters = poseAndShapeParameters.rotationParameters

      // assuming all parameters are independent, we can sum the log probabilities
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

  //--------------------------------------------------------------------------------------------
  // Likelihood
  //--------------------------------------------------------------------------------------------

  /** Likelihood model that the closest point of the target mesh is observed up
    * to Gaussian noise.
    */
  class ClosestPointEvaluator(
      fullModel: PointDistributionModel[_3D, TriangleMesh],
      rotationCenter: Point[_3D],
      targetMesh: TriangleMesh[_3D]
  )(implicit val rng: scalismo.utils.Random)
      extends MHDistributionEvaluator[Parameters] {

    // decimation of the model is needed for speedup
    val model = fullModel.newReference(
      fullModel.reference.operations.decimate(250),
      TriangleMeshInterpolator3D()
    )
    override def logValue(sample: MHSample[Parameters]): Double = {
      val noiseModel =
        breeze.stats.distributions.MultivariateGaussian(
          DenseVector.zeros[Double](3),
          DenseMatrix.eye[Double](3) * Math.abs(sample.parameters.noiseStddev * sample.parameters.noiseStddev)
        )

      val poseAndShapeParameters = sample.parameters.poseAndShapeParameters
      val translationParameters = poseAndShapeParameters.translationParameters
      val rotationParameters = poseAndShapeParameters.rotationParameters

      val poseTransformation = Parameters.poseTransformationForParameters(
        translationParameters,
        rotationParameters,
        rotationCenter
      )

      val transformedModel = model
        .instance(poseAndShapeParameters.shapeParameters.coefficients)
        .transform(poseTransformation)

      val logValues = for (pointOnTarget <- targetMesh.pointSet.points) yield {
        val u =
          transformedModel.operations
            .closestPointOnSurface(pointOnTarget)
            .point - pointOnTarget
        noiseModel.logPdf(u.toBreezeVector)
      }

      logValues.sum
    }
  }


  // A simple evaluator that accepts all samples - for debugging purpose only
  //
  case class AcceptAllEvaluator() extends MHDistributionEvaluator[Parameters] {
    override def logValue(sample: MHSample[Parameters]): Double = {
      0.0
    }
  }


  // --------------------------------------------------------------------------------------------
  //  Main fitting method
  // --------------------------------------------------------------------------------------------
  def fitShape(
      model: PointDistributionModel[_3D, TriangleMesh],
      rotationCenter: Point[_3D],
      targetMesh: TriangleMesh[_3D],
      initialParameters: PoseAndShapeParameters,
      shapeModelTransformationView: ShapeModelTransformationView
  )(implicit  rng: scalismo.utils.Random ): Unit = {

    //
    // Setting up the target distribtion
    //
    val prior = PriorEvaluator(model)
    val likelihood =
      new ClosestPointEvaluator(model, rotationCenter, targetMesh).cached
    val posterior = ProductEvaluator(prior, likelihood)

    //
    // Define the proposals
    //
    val rotationProposal = MHProductProposal(
      GaussianRandomWalkProposal(0.01, "rho").forType[Double],
      GaussianRandomWalkProposal(0.01, "phi").forType[Double],
      GaussianRandomWalkProposal(0.01, "psi").forType[Double]
    ).forType[RotationParameters]

    val translationProposal = MHProductProposal(
      GaussianRandomWalkProposal(0.1, "tx").forType[Double],
      GaussianRandomWalkProposal(0.1, "ty").forType[Double],
      GaussianRandomWalkProposal(0.1, "tz").forType[Double]
    ).forType[TranslationParameters]

    val shapeProposalLeading =
      GaussianRandomWalkProposal(0.05, "shape-0-5")
        .partial(0 until 5)
        .forType[ShapeParameters]
    val shapeProposalRemaining =
      GaussianRandomWalkProposal(0.05, "shape-6-")
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
      ).forType[PoseAndShapeParameters]
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
        .relabel("shape-only-leading")

    val poseAndShapeRemainingShapeOnlyProposal =
      MHProductProposal(
        identTranslationProposal,
        identRotationProposal,
        shapeProposalRemaining
      )
        .forType[PoseAndShapeParameters]
        .relabel("shape-only-trailing")

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
        .relabel("noise-only")

    val poseAndShapeOnlyProposal =
      MHProductProposal(mixturePoseAndShapeProposal, identNoiseProposal)
        .forType[Parameters]

    val fullproposal = MHMixtureProposal(
      (0.9, poseAndShapeOnlyProposal),
      (0.0, noiseOnlyProposal)
    )

    //
    // setting up the sampler
    //
    val logger = MHSampleLogger[Parameters]()
    val chain = MetropolisHastings(fullproposal, posterior)

    val initialParameters = Parameters(
      PoseAndShapeParameters(
        TranslationParameters(EuclideanVector3D(0, 0, 0)),
        RotationParameters((0.0, 0.0, 0.0)),
        ShapeParameters(DenseVector.zeros[Double](model.rank))
      ),
      noiseStddev = 5.0
    )

    val mhIterator =
      chain.iterator(MHSample(initialParameters, "inital"), logger)

    //
    // visualizing the state of the iteration
    //
    val samplingIterator =
      for ((sample, iteration) <- mhIterator.zipWithIndex) yield {
        println("iteration " + iteration)
        if (iteration % 500 == 0) {
          val poseAndShapeParameters = sample.parameters.poseAndShapeParameters
          val poseTransformation = Parameters.poseTransformationForParameters(
            poseAndShapeParameters.translationParameters,
            poseAndShapeParameters.rotationParameters,
            rotationCenter
          )
          shapeModelTransformationView.shapeTransformationView.coefficients =
            poseAndShapeParameters.shapeParameters.coefficients
          shapeModelTransformationView.poseTransformationView.transformation =
            poseTransformation

        }
        sample
      }

    // obtaining the samples
    val samples = samplingIterator.drop(1000).take(5000).toIndexedSeq

    // extracting and visualizing the most likely solution
    val bestSample = samples.maxBy(posterior.logValue)
    val sampleGroup = ui.createGroup("sampleGroup")

    val bestSurfaceReconstruction = model
      .instance(
        bestSample.parameters.poseAndShapeParameters.shapeParameters.coefficients
      )
      .transform(
        Parameters.poseTransformationForParameters(
          bestSample.parameters.poseAndShapeParameters.translationParameters,
          bestSample.parameters.poseAndShapeParameters.rotationParameters,
          rotationCenter
        )
      )
    ui.show(sampleGroup, bestSurfaceReconstruction, "best sample")

    // Some diagnostics and statistics

    println(
      "acceptance ratio " + logger.samples.takeLast(1000).acceptanceRatios
    )

    val dataframe = DataFrame.fromColumns(Seq(
      Column.ofContinuous(samples.map(s => posterior.logValue(s)), "logPosterior"),      
      Column.ofContinuous(samples.map(s => s.parameters.poseAndShapeParameters.shapeParameters.coefficients(0)), "shape-0"),
      Column.ofContinuous(samples.map(s => s.parameters.poseAndShapeParameters.translationParameters.translationVector.x), "tx"),
      Column.ofContinuous(samples.map(s => s.parameters.poseAndShapeParameters.translationParameters.translationVector.y), "ty"),
      Column.ofContinuous(samples.map(s => s.parameters.poseAndShapeParameters.translationParameters.translationVector.z), "tz"),
      Column.ofContinuous(samples.map(s => s.parameters.noiseStddev), "noise")
    ))

    dataframe.plot.tracePlot("logPosterior", "Trace plot").show()
    dataframe.plot.pairPlot(Seq("shape-0", "tx", "ty", "tz", "noise"), "pairs").show()    



  }

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val model = StatisticalModelIO
      .readStatisticalTriangleMeshModel3D(
        new java.io.File("datasets/femur/ssm.h5")
      )
      .get

    // val model = StatisticalModelIO
    //   .readStatisticalTriangleMeshModel3D(new java.io.File("datasets/femur/gpmodel.h5"))
    //   .get

    // fake data
    val fakePoseTransform = Parameters.poseTransformationForParameters(
      TranslationParameters(EuclideanVector3D(30, 0, 10)),
      RotationParameters((0.1, 0.3, -0.1)),
      Point(0, 0, 0)
    )

    val fakeCoeffs = DenseVector.zeros[Double](model.rank)
    fakeCoeffs(0) = 3.0
    fakeCoeffs(1) = 2.0
    fakeCoeffs(2) = 1.0

    // val targetMesh = model.instance(fakeCoeffs).transform(fakePoseTransform)

    val targetMesh =
      MeshIO
        .readMesh(new java.io.File("datasets/femur/fragments/fragment-9.stl"))
        .get

    val modelGroup = ui.createGroup("model")
    val modelView = ui.show(modelGroup, model, "model")

    val targetGroup = ui.createGroup("target")
    ui.show(targetGroup, targetMesh, "target")

    def centerOfMass(mesh: TriangleMesh[_3D]): Point[_3D] = {
      val normFactor = 1.0 / mesh.pointSet.numberOfPoints
      val origin = Point.origin[_3D]
      mesh.pointSet.points.foldLeft(origin)((sum, point) =>
        sum + (point - origin) * normFactor
      )
    }

    val rotationCenter = model.mean.pointSet.centerOfMass

    val initialParameters = PoseAndShapeParameters(
      TranslationParameters(EuclideanVector3D(0, 0, 0)),
      RotationParameters(0, 0, 0),
      ShapeParameters(DenseVector.zeros[Double](model.rank))
    )

    fitShape(
      model,
      rotationCenter,
      targetMesh,
      initialParameters,
      modelView.shapeModelTransformationView
    )

  }
}
