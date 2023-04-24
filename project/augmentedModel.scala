//> using scala "3.1.2"
//> using lib "ch.unibas.cs.gravis::scalismo-ui:0.91.0"
import scalismo.geometry._
import scalismo.common._
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.mesh._
import scalismo.io.{StatisticalModelIO, MeshIO}
import scalismo.statisticalmodel._
import scalismo.numerics.UniformMeshSampler3D
import scalismo.kernels._

import scalismo.ui.api._

import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.io.LandmarkIO

object augmentedModel extends App {

    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    // Create Directory to store results in
    val resultDir = new java.io.File("results/lowRankSamples")
    resultDir.mkdirs()

    // load reference-mesh
    val referenceMesh = MeshIO
        .readMesh(
            new java.io.File("project-data/reference-mesh/reference.stl")
        )
        .get

    // Load fitted model
    val pcaModel = StatisticalModelIO.readStatisticalTriangleMeshModel3D(new java.io.File("fittedModel.h5")).get
    val gpSSM = pcaModel.gp.interpolate(TriangleMeshInterpolator3D())

    val covSSM : MatrixValuedPDKernel[_3D] = gpSSM.cov

    val scalarValuedGaussianKernel1 : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0, scaleFactor = 0.5)
    val scalarValuedGaussianKernel2 : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0, scaleFactor = 0.5)
    val scalarValuedGaussianKernel3 : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0, scaleFactor = 1)

    val matrixValuedGaussianKernel1 = DiagonalKernel3D(scalarValuedGaussianKernel1, scalarValuedGaussianKernel2, scalarValuedGaussianKernel3)

    val matrixValuedGaussianKernel  = matrixValuedGaussianKernel1 

    val augmentedCov = covSSM + matrixValuedGaussianKernel

    val augmentedGP = GaussianProcess(gpSSM.mean, augmentedCov)

    val lowRankAugmentedGP = LowRankGaussianProcess.approximateGPCholesky(
    referenceMesh,
    augmentedGP,
    relativeTolerance = 0.01,
    interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]()
    )
    val augmentedSSM = PointDistributionModel3D(pcaModel.reference, lowRankAugmentedGP)
    
    val refGroup = ui.createGroup("reference")
    val testview = ui.show(refGroup, augmentedSSM,"test")
    val ripview = ui.show(refGroup, augmentedSSM.reference,"testref")
    // Save model
    val augmentedModel = StatisticalModelIO.writeStatisticalTriangleMeshModel3D(augmentedSSM,
        new java.io.File("augmentedModel.h5")) // saving the model
    printf("done")
}
