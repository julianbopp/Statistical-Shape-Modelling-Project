//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.92-SNAPSHOT"
//> using lib "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.io.LandmarkIO
import scalismo.io.MeshIO
import scalismo.geometry._3D
import scalismo.registration.LandmarkRegistration
import scalismo.geometry.{Point3D, EuclideanVector3D}
import scalismo.ui.api.ScalismoUI
import scalismo.geometry.{Landmark, Landmark3D}

/** 
 * Implementation of generalized Procrustes analysis
  */
object GPA {

  implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42L)

  case class MeanAndAlignedLandmarks(
      meanLandmarks: Seq[Landmark[_3D]],
      alignedLandmarks: Seq[Seq[Landmark[_3D]]]
  )

  def doGPA(
      targetLandmarks: Seq[Landmark[_3D]],
      landmarks: Seq[Seq[Landmark[_3D]]],
      numberOfIterations: Int
  ): MeanAndAlignedLandmarks = {

    if (numberOfIterations > 0) then {
      val meanLandmarks = computeMeanLandmarks(landmarks)
      val alignedLandmarks: Seq[Seq[Landmark[_3D]]] =
        alignLandmarks(meanLandmarks, landmarks)
      doGPA(meanLandmarks, alignedLandmarks, numberOfIterations - 1)
    } else {
      MeanAndAlignedLandmarks(targetLandmarks, landmarks)
    }
  }

  def alignLandmarks(
      target: Seq[Landmark[_3D]],
      landmarkSeq: Seq[Seq[Landmark[_3D]]]
  ): Seq[Seq[Landmark[_3D]]] = {
    for (landmarks <- landmarkSeq) yield {
      val transformation = LandmarkRegistration.rigid3DLandmarkRegistration(
        landmarks,
        target,
        Point3D(0, 0, 0)
      )
      landmarks.map(landmark => landmark.transform(transformation))
    }
  }

  def computeLandmarkMean(landmarks: Seq[Landmark[_3D]]): Landmark[_3D] = {
    
    var avgVec = EuclideanVector3D(0, 0, 0)
    for (landmark <- landmarks) yield {
        avgVec += landmark.point - Point3D(0, 0, 0)
    }
    avgVec = avgVec / landmarks.length
    Landmark3D(landmarks.head.id, avgVec.toPoint)
  }

  def computeMeanLandmarks(
      landmarkSeq: Seq[Seq[Landmark[_3D]]]
  ): Seq[Landmark[_3D]] = {

    // the following code extracts each corresponding landmark and computes the mean.
    // The code could, of course, be rewritten more succinctly using a for loop

    val landmark0 = computeLandmarkMean(
      landmarkSeq.map(landmarks => landmarks(0))
    )
    val landmark1 = computeLandmarkMean(
      landmarkSeq.map(landmarks => landmarks(1))
    )
    val landmark2 = computeLandmarkMean(
      landmarkSeq.map(landmarks => landmarks(2))
    )
    val landmark3 = computeLandmarkMean(
      landmarkSeq.map(landmarks => landmarks(3))
    )
    val landmark4 = computeLandmarkMean(
      landmarkSeq.map(landmarks => landmarks(4))
    )
    val landmark5 = computeLandmarkMean(
      landmarkSeq.map(landmarks => landmarks(5))
    )

    Seq(landmark0, landmark1, landmark2, landmark3, landmark4, landmark5)
  }

  def main(args: Array[String]): Unit = {

    val dataDir = new java.io.File("project-data")

    val ui = ScalismoUI()

    val refGroup = ui.createGroup("reference")
    val alignedGroup = ui.createGroup("mean-landmarks")

    val referenceLandmarks = LandmarkIO
      .readLandmarksJson3D(
        new java.io.File(dataDir, "/reference-landmarks/reference.json")
      )
      .get

    for landmark <- referenceLandmarks do
      ui.show(refGroup, landmark, landmark.id)

    val fileIds = 0 until 46

    // loop through all the ids
    val landmarksSeq = for (fileId <- fileIds) yield {
      val landmarkFile = new java.io.File(dataDir, s"landmarks/$fileId.json")
      LandmarkIO.readLandmarksJson3D(landmarkFile).get
    }

    val MeanAndAlignedLandmarks(meanLandmarks, alignedLandmarks) =
      doGPA(referenceLandmarks, landmarksSeq, 10)

    println(meanLandmarks)

    for landmark <- meanLandmarks do
      
      ui.show(alignedGroup, landmark, landmark.id)

  }

}
