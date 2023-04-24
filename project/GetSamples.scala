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

object GetSamples extends App {

    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    // Create Directory to store results in
    val resultDir = new java.io.File("results/fittedSamples")
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

    
    val fittedmodel = StatisticalModelIO.readStatisticalTriangleMeshModel3D(new java.io.File("fittedModel.h5")).get
    val gp = fittedmodel.gp

    val interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]()
    val contGP = fittedmodel.gp.interpolate(interpolator)

    val lowRankGP = LowRankGaussianProcess.approximateGPCholesky(referenceMesh,contGP,relativeTolerance = 0.07,interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]())

    val numOfSamples: Int = 46 
    val sampleGroup = ui.createGroup("gp-sample")

    // Sample from the LowRankGaussianProcess object
    val defFieldSamples = (0 until numOfSamples).map(_ => lowRankGP.sample()).toList
    //val landmarkViews = ui.filter[LandmarkView](group, (v : LandmarkView) => true)

    for (i <- 0 until numOfSamples) {
        val defField = defFieldSamples(i) 
        val sampleMesh = referenceMesh.copy().transform((p : Point[_3D]) => p + defField(p))

        val p0 = referenceMesh.pointSet.point(PointId(961))
        val p1 = referenceMesh.pointSet.point(PointId(1307))
        val p2 = referenceMesh.pointSet.point(PointId(1526))
        val p3 = referenceMesh.pointSet.point(PointId(3333))
        val p4 = referenceMesh.pointSet.point(PointId(3785))
        val p5 = referenceMesh.pointSet.point(PointId(2509))

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
    print(model.mean)
    ui.show(modelGroup , model , "gp model")
    printf("Done")
}
