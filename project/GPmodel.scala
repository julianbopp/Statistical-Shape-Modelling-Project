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

object GPmodel extends App {

    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    // Create Directory to store results in
    val resultDir = new java.io.File("results/lowRankSamples")
    resultDir.mkdirs()

    // Load reference mesh and display in separate group
    val referenceMesh = MeshIO.readMesh(new java.io.File("project-data/reference-mesh/reference.stl")).get
    val modelGroup = ui.createGroup("gp-model")
    val referenceView = ui.show(modelGroup, referenceMesh, "reference")

    // Define mean deformation as zero mean (i.e. we assume reference shape to be mean shape)
    val zeroMean = Field(EuclideanSpace3D, (pt:Point[_3D]) => EuclideanVector3D(0,0,0))

    // Define Gaussian Kernel
    val scalarValuedGaussianKernel : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0)
    val matrixValuedGaussianKernel = DiagonalKernel3D(scalarValuedGaussianKernel, 3) 

    // Define Gaussian process
    val gp = GaussianProcess3D[EuclideanVector[_3D]](zeroMean, matrixValuedGaussianKernel)

   // Define Low-rank approximation 
    val lowRankGP = LowRankGaussianProcess.approximateGPCholesky(
        referenceMesh,
        gp,
        relativeTolerance = 0.1,
        interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]()
    )

    var  defField : Field[_3D, EuclideanVector[_3D]]= lowRankGP.sample()
    referenceMesh.transform((p : Point[_3D]) => p + defField(p))
    var pdm = PointDistributionModel3D(referenceMesh, lowRankGP)
    var pdmView = ui.show(modelGroup, pdm, "group")


    val numOfSamples: Int = 20
    val sampleGroup = ui.createGroup("gp-sample")

    // Sample from the LowRankGaussianProcess object
    val defFieldSamples = (0 until numOfSamples).map(_ => lowRankGP.sample()).toList
    //val landmarkViews = ui.filter[LandmarkView](group, (v : LandmarkView) => true)

    for (i <- 0 until numOfSamples) {
        val defField = defFieldSamples(i) 
        val sampleMesh = referenceMesh.transform((p : Point[_3D]) => p + defField(p))
        ui.show(sampleGroup, sampleMesh, s"sample-$i")
        MeshIO
            .writeMesh(sampleMesh, new java.io.File(resultDir, s"$i.stl"))
            .get
    }

}
