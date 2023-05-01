//> using scala "3.1.2"
//> using lib "ch.unibas.cs.gravis::scalismo-ui:0.91.0"
import scalismo.ui.api._

import scalismo.geometry._
import scalismo.common._
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.mesh._
import scalismo.io.{StatisticalModelIO, MeshIO, LandmarkIO}
import scalismo.statisticalmodel._
import scalismo.registration._
import scalismo.statisticalmodel.dataset._
import scalismo.numerics.PivotedCholesky.RelativeTolerance

object ICP_model extends App {
    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    val dsGroup = ui.createGroup("datasets")

    val meshFiles = new java.io.File("results/fitted/").listFiles
    val (meshes, meshViews) = meshFiles
    .map(meshFile => {
        val mesh = MeshIO.readMesh(meshFile).get
        val meshView = ui.show(dsGroup, mesh, "mesh")
        (mesh, meshView) // return a tuple of the mesh and the associated view
    })
    .unzip // take the tuples apart, to get a sequence of meshes and one of meshViews

    // Load the reference mesh
    val reference = MeshIO.readMesh(new java.io.File("project-data/reference-mesh/reference.stl")).get
    val toAlign: IndexedSeq[TriangleMesh[_3D]] = meshes

    // Load the Landmarks, these pointIds correspond to the Landmarks L0,...,L5
    val pointIds = IndexedSeq(961, 1307, 1526, 3333, 3785, 2509)
    val refLandmarks = pointIds.map { id =>
    Landmark(s"L_$id", reference.pointSet.point(PointId(id)))
    }
    // Save the reference
    LandmarkIO.writeLandmarksJson[_3D](refLandmarks, new java.io.File("results/fitted-reference-landmarks/reference.json")).get

    val alignedMeshes = toAlign.map { mesh =>
    val landmarks = pointIds.map { id => Landmark("L_" + id, mesh.pointSet.point(PointId(id))) }
    val rigidTrans = LandmarkRegistration.rigid3DLandmarkRegistration(
        landmarks,
        refLandmarks,
        center = Point(0, 0, 0)
    )
    mesh.transform(rigidTrans)
    }

    val dc = DataCollection.fromTriangleMesh3DSequence(reference, meshes)

    // Create PDM model using PCA
    val modelFromDataCollection = PointDistributionModel.createUsingPCA(dc)

    // Save model
    val fittedModel = StatisticalModelIO.writeStatisticalTriangleMeshModel3D(modelFromDataCollection,
        new java.io.File("fittedModel.h5")) // saving the model

    val modelGroup = ui.createGroup("modelGroup")
    ui.show(modelGroup, modelFromDataCollection, "ModelDC")
    


}
