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
import scalismo.plot.data.DataFrame
import scalismo.plot.data.DataFrame.Column
import scalismo.plot.plottarget.PlotTargets.plotTargetBrowser

/** 
 * Performs length and width measurements on a set of femur shapes. The length and with is 
 * determined by the corresponding landmarks
 */
object LengthMeasurements {

  case class Measurement(id: Int, length: Double, width: Double)

  case class MeanAndVariance(mean: Double, variance: Double)

  def meanAndVariance(measurements: Seq[Double]): MeanAndVariance = {
    val mean = measurements.sum / measurements.length

    val variance = measurements.map(measurement => (measurement - mean) * (measurement - mean)).sum * 1.0 / measurements.length
    
    MeanAndVariance(mean, variance)
  }

  def main(args: Array[String]): Unit = {
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42L)

    val dataDir = new java.io.File("project-data")

    val referenceLandmarks = LandmarkIO
      .readLandmarksJson3D(
        new java.io.File(dataDir, "/reference-landmarks/reference.json")
      )
      .get

    val fileIds = 0 until 46

    // Computing the measurements
    val measurements = for (fileId <- fileIds) yield {
      println(s"processing $fileId")

      val landmarkFile =
        new java.io.File(dataDir, s"landmarks/$fileId.json")

      val landmarks = LandmarkIO.readLandmarksJson3D(landmarkFile).get
      val lmL2 = landmarks.find(lm => lm.id == "L2").get
      val lmL5 = landmarks.find(lm => lm.id == "L5").get

      val lmL3 = landmarks.find(lm => lm.id == "L3").get
      val lmL4 = landmarks.find(lm => lm.id == "L4").get

      val length = (lmL5.point - lmL2.point).norm
      val width = (lmL3.point - lmL4.point).norm

      Measurement(fileId, length, width)
    }

    println("mean and variance length: " + meanAndVariance(measurements.map(m => m.length)))
    println("mean and variance width: " + meanAndVariance(measurements.map(m => m.width)))

    val dataFrame = DataFrame(
      Seq(
        Column.ofContinuous(measurements.map(m => m.length), "length"),
        Column.ofContinuous(measurements.map(m => m.width), "width")
      )
    )

    dataFrame.plot.scatterPlot("length", "width", "Length-vs-Width").show()
    lowRankMeasurements()
  }

  def lowRankMeasurements(): Unit = {

    val referenceLandmarks = LandmarkIO
      .readLandmarksJson3D(
        new java.io.File("project-data/reference-landmarks/reference.json")
      )
      .get

    val fileIds = 0 until 20

    // Computing the measurements
    val measurements = for (fileId <- fileIds) yield {
      println(s"processing $fileId")

      val landmarkFile =
        new java.io.File(s"results/lowRankSamples/$fileId.json")

      val landmarks = LandmarkIO.readLandmarksJson3D(landmarkFile).get
      val lmL2 = landmarks.find(lm => lm.id == "L2").get
      val lmL5 = landmarks.find(lm => lm.id == "L5").get

      val lmL3 = landmarks.find(lm => lm.id == "L3").get
      val lmL4 = landmarks.find(lm => lm.id == "L4").get

      val length = (lmL5.point - lmL2.point).norm
      val width = (lmL3.point - lmL4.point).norm

      Measurement(fileId, length, width)
    }

    println("mean and variance length: " + meanAndVariance(measurements.map(m => m.length)))
    println("mean and variance width: " + meanAndVariance(measurements.map(m => m.width)))

    val dataFrame = DataFrame(
      Seq(
        Column.ofContinuous(measurements.map(m => m.length), "length"),
        Column.ofContinuous(measurements.map(m => m.width), "width")
      )
    )

    dataFrame.plot.scatterPlot("length", "width", "Length-vs-Width").show()
  }
  
}
