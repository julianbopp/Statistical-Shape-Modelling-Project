//> using scala "3.2"
//> using repository "sonatype:snapshots"
////> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using lib "ch.unibas.cs.gravis::scalismo-ui:0.91.0"
//> using lib "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.io.LandmarkIO
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.geometry._
import scalismo.common._
import scalismo.common.interpolation.{NearestNeighborInterpolator, TriangleMeshInterpolator3D}
import scalismo.registration.LandmarkRegistration
import scalismo.ui.api._
import scalismo.transformations._
import scalismo.statisticalmodel._
import scalismo.kernels._
import java.awt.Color
import java.io.File
import scala.swing.Color

import breeze.linalg.DenseVector



object GPmodel extends App{
    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    val dataDir = new java.io.File("project-data")

    // Load reference mesh
    val referenceMesh = MeshIO.readMesh(new java.io.File(dataDir, "reference-mesh/reference.stl")).get

    // Load reference landmarks
    val referenceLandmarks = LandmarkIO.readLandmarksJson3D(new java.io.File(dataDir, "/reference-landmarks/reference.json")).get

    val modelGroup = ui.createGroup("gp-model")
    val refView  = ui.show(modelGroup, referenceMesh, "reference")
    val rigidFemurs = ui.createGroup("RigidFemurs")
    val alignedFemurs = ui.createGroup("RigidFemurs")

    //Initialize mean for Gaussian Process
    val zeroMean = Field(EuclideanSpace3D, (pt: Point[_3D]) => EuclideanVector3D(0, 0, 0))

    // Create result directory
    val resultDir = new java.io.File("results/aligned")
    resultDir.mkdirs() // dir does not exist create it.

    //Load a rigid mesh
    val meshRigid = MeshIO.readMesh(new java.io.File(dataDir, "meshes/30.stl")).get 

    //Load a rigid landmark
    val refLandmarks = LandmarkIO.readLandmarksJson3D(new java.io.File(dataDir, "/landmarks/30.json")).get

    // Rigid alignment

    val bestTransform: RigidTransformation[ _3D ] = LandmarkRegistration.rigid3DLandmarkRegistration(refLandmarks, referenceLandmarks, center = Point(0, 0, 0))
    val alignedFemur = meshRigid.transform(bestTransform)
    val alignedLandmarks = refLandmarks.map(lm => lm.copy(point = bestTransform(lm.point)))
    val rigidFemurView = ui.show(rigidFemurs, meshRigid, "Rigid_Femur")
    val alignedFemurView = ui.show(alignedFemurs, alignedFemur, "Aligned_Femur")

    rigidFemurView.color = java.awt.Color.red
    alignedFemurView.color = java.awt.Color.green

    val refMeshColor = new Color(115,238,70)
    val defMeshColor = new Color(101,70,108)
    



    //setup kernels
    val scalarValuedGaussianKernel1 : PDKernel[_3D]= GaussianKernel3D(sigma = 40.0)
    val scalarValuedGaussianKernel2 : PDKernel[_3D]= GaussianKernel3D(sigma = 60.0)
    val scalarValuedGaussianKernel3 : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0)

    val matrixValuedGaussianKernel1 = DiagonalKernel3D(scalarValuedGaussianKernel1, scalarValuedGaussianKernel2, scalarValuedGaussianKernel3)

    val scalarValuedGaussianKernel4 : PDKernel[_3D]= GaussianKernel3D(sigma = 50.0)
    val scalarValuedGaussianKernel5 : PDKernel[_3D]= GaussianKernel3D(sigma = 60.0)
    val scalarValuedGaussianKernel6 : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0)

    val matrixValuedGaussianKernel2 = DiagonalKernel3D(scalarValuedGaussianKernel4, scalarValuedGaussianKernel5, scalarValuedGaussianKernel6)

    val matrixValuedGaussianKernel  = matrixValuedGaussianKernel1 + matrixValuedGaussianKernel2

    val gp = GaussianProcess3D[EuclideanVector[_3D]](zeroMean, matrixValuedGaussianKernel)

    val lowRankGP = LowRankGaussianProcess.approximateGPCholesky(referenceMesh,gp,relativeTolerance = 0.07,interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]())

    val deformedGroup = ui.createGroup("deformed")
    val refMeshView2  = ui.show(deformedGroup, referenceMesh, "referenceMesh")


    //refMeshView2.color = defMeshColor
    val gpDefView = ui.addTransformation(deformedGroup, lowRankGP, "RefMeshDeformedByGp")

    val model = PointDistributionModel(referenceMesh,lowRankGP)  //creating the model
     ui.show(modelGroup , model , "gp model")

    val gpModel = StatisticalModelIO.writeStatisticalTriangleMeshModel3D(model, new java.io.File(dataDir, "/gpmodel.h5")) //saving the model







}