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

    // Load reference Landmarks
    val referenceLandmarks = LandmarkIO
      .readLandmarksJson3D(
        new java.io.File("project-data/reference-landmarks/reference.json")
      )
      .get

    // Define mean deformation as zero mean (i.e. we assume reference shape to be mean shape)
    val zeroMean = Field(EuclideanSpace3D, (pt:Point[_3D]) => EuclideanVector3D(0,0,0))

    // Setup kernels
    val scalarValuedGaussianKernel1 : PDKernel[_3D]= GaussianKernel3D(sigma = 50.0, scaleFactor = 5)
    val scalarValuedGaussianKernel2 : PDKernel[_3D]= GaussianKernel3D(sigma = 50.0, scaleFactor = 5)
    val scalarValuedGaussianKernel3 : PDKernel[_3D]= GaussianKernel3D(sigma = 200.0, scaleFactor = 20)

    val matrixValuedGaussianKernel1 = DiagonalKernel3D(scalarValuedGaussianKernel1, scalarValuedGaussianKernel2, scalarValuedGaussianKernel3)

    val matrixValuedGaussianKernel  = matrixValuedGaussianKernel1 

    val gp = GaussianProcess3D[EuclideanVector[_3D]](zeroMean, matrixValuedGaussianKernel)

    val lowRankGP = LowRankGaussianProcess.approximateGPCholesky(referenceMesh,gp,relativeTolerance = 0.07,interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]())


    val numOfSamples: Int = 47 
    val sampleGroup = ui.createGroup("gp-sample")

    // Sample from the LowRankGaussianProcess object
    val defFieldSamples = (0 until numOfSamples).map(_ => lowRankGP.sample()).toList
    //val landmarkViews = ui.filter[LandmarkView](group, (v : LandmarkView) => true)

    for (i <- 0 until numOfSamples) {
        val defField = defFieldSamples(i) 
        val sampleMesh = referenceMesh.copy().transform((p : Point[_3D]) => p + defField(p))

        val p0 = referenceLandmarks(0).point
        val p1 = referenceLandmarks(1).point
        val p2 = referenceLandmarks(2).point
        val p3 = referenceLandmarks(3).point
        val p4 = referenceLandmarks(4).point
        val p5 = referenceLandmarks(5).point

        val tp0 = p0 + defField(p0)
        val tp1 = p1 + defField(p1)
        val tp2 = p2 + defField(p2)
        val tp3 = p3 + defField(p3)
        val tp4 = p4 + defField(p4)
        val tp5 = p5 + defField(p5)

        val L0 = Landmark3D("L0", tp0)
        val L1 = Landmark3D("L1", tp1)
        val L2 = Landmark3D("L2", tp2)
        val L3 = Landmark3D("L3", tp3)
        val L4 = Landmark3D("L4", tp4)
        val L5 = Landmark3D("L5", tp5)
        
        val sampleLandmarks = Seq(L0,L1,L2,L3,L4,L5)

        ui.show(sampleGroup, sampleMesh, s"sample-$i")
        MeshIO
            .writeMesh(sampleMesh, new java.io.File(resultDir, s"$i.stl"))
            .get
        LandmarkIO
        .writeLandmarksJson[_3D](
            sampleLandmarks,
            new java.io.File(resultDir, s"$i.json")
        )
        .get
    }
    val model = PointDistributionModel(referenceMesh,lowRankGP)  //creating the model
    ui.show(modelGroup , model , "gp model")

    val gpModel = StatisticalModelIO.writeStatisticalTriangleMeshModel3D(model, new java.io.File("gpmodel.h5")) //saving the model


    val p1 = referenceLandmarks(3).point
    val p2 = referenceLandmarks(4).point
    print((p1-p2).norm)
}
