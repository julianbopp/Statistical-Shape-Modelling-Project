//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using lib "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.io.LandmarkIO
import scalismo.io.MeshIO
import scalismo.geometry._3D
import scalismo.registration.LandmarkRegistration
import scalismo.geometry.Point3D
import scalismo.ui.api.ScalismoUI
import scalismo.geometry.Landmark


/** Rigid alignment of the data as required in exercise sheet 1
  */
object RigidAlignment extends App {

  implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42L)
  val ui = ScalismoUI()


  val dataDir = new java.io.File("project-data")

  // load reference mesh
  val referenceMesh = MeshIO
    .readMesh(
      new java.io.File(dataDir, "reference-mesh/reference.stl")
    )
    .get
  val referenceLandmarks = LandmarkIO
    .readLandmarksJson3D(
      new java.io.File(dataDir, "/reference-landmarks/reference.json")
    )
    .get


  val refGroup = ui.createGroup("reference")
  val refView = ui.show(refGroup, referenceMesh, "reference")  

  println("the center of mass of the reference shape is " + referenceMesh.pointSet.centerOfMass)

  // create result directory
  val resultDir = new java.io.File("results/aligned")
  resultDir.mkdirs()

  val fileIds = 0 until 46 // change to 46 to process all the files

  val alignedGroup = ui.createGroup("aligned meshes")

  // loop through all the ids and load the files
  for (fileId <- fileIds) {
    println(s"processing $fileId")
    
    // loading data for id
    val landmarkFile =
      new java.io.File(dataDir, s"landmarks/$fileId.json")
    
    val meshFile =
      new java.io.File(dataDir, s"meshes/$fileId.stl")
    val landmarks = LandmarkIO.readLandmarksJson3D(landmarkFile).get
    val mesh = MeshIO.readMesh(meshFile).get

    // compute transformation transform mesh and landmarks
    val transformation = LandmarkRegistration.rigid3DLandmarkRegistration(
      landmarks,
      referenceLandmarks,
      Point3D(0, 0, 0)
    )
    
    // apply transformation to meshes and landmarks
    val transformedMesh = mesh.transform(transformation)
    val transformedLandmarks = landmarks.map(_.transform(transformation))

    // Write everything into the results directory
    MeshIO
      .writeMesh(transformedMesh, new java.io.File(resultDir, s"$fileId.stl"))
      .get
    LandmarkIO
      .writeLandmarksJson[_3D](
        transformedLandmarks,
        new java.io.File(resultDir, s"$fileId.json")
      )
      .get

    ui.show(alignedGroup, transformedMesh, "transformed")
  }
}