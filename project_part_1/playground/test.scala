//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using lib "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"


// Basic geometric primitives
import scalismo.geometry.{_3D, Point, Point3D}
import scalismo.geometry.{EuclideanVector}
import scalismo.geometry.{IntVector, IntVector3D} 
import scalismo.geometry.Landmark

import scalismo.common.PointId

// Geometric objects
import scalismo.mesh.TriangleMesh
import scalismo.mesh.TriangleId
import scalismo.image.{DiscreteImage, DiscreteImage3D}
import scalismo.statisticalmodel.PointDistributionModel 

// IO Methods
import scalismo.io.ImageIO; 
import scalismo.io.StatisticalModelIO
import scalismo.io.{MeshIO, StatisticalModelIO}

// Visualization
import scalismo.ui.api.ScalismoUI
import scalismo.ui.api.LandmarkView
import scalismo.io.LandmarkIO

object test extends App {

    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    val mesh = MeshIO.readMesh(new java.io.File("/home/bobby/git/Statistical-Shape-Modelling-Project/project/project-data/reference-mesh/reference.stl")).get
    val landmarks = LandmarkIO.readLandmarksJson3D(new java.io.File("/home/bobby/git/Statistical-Shape-Modelling-Project/project/project-data/reference-landmarks/reference.json")).get
    val fittedmodel = StatisticalModelIO.readStatisticalTriangleMeshModel3D(new java.io.File("../project/fittedModel.h5")).get

    val group = ui.createGroup("testGroup")
    val modelgroup = ui.createGroup("modelGroup")
    val meshView = ui.show(group, mesh, "testView")
    ui.show(group, landmarks, "landmarks")
    ui.show(modelgroup, fittedmodel, "model")

    val meshtest = MeshIO.readMesh(new java.io.File("/home/bobby/git/Statistical-Shape-Modelling-Project/project/results/fittedSamples/0.stl")).get
    val landmarkstest = LandmarkIO.readLandmarksJson3D(new java.io.File("/home/bobby/git/Statistical-Shape-Modelling-Project/project/results/fittedSamples/0.json")).get
    ui.show(group, meshtest, "testmesh")
    ui.show(group, landmarkstest, "testllandmarks")
}
