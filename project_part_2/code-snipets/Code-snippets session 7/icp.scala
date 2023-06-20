//> using scala "3.2.1"
//> using repository "sonatype:snapshots"
//> using lib "ch.unibas.cs.gravis::scalismo-ui:0.92-RC1"
//> using lib "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.kernels.*
import scalismo.statisticalmodel.{LowRankGaussianProcess, GaussianProcess}
import scalismo.common.*
import scalismo.geometry.*
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.io.{MeshIO, LandmarkIO}

import scalismo.ui.api.ScalismoUI

import scalismo.plot.data.*
import scalismo.plot.plottarget.PlotTargets.plotTargetBrowser
import scalismo.plot.data.DataFrame.Column

import breeze.stats.distributions.Rand.FixedSeed.randBasis
import scalismo.common.interpolation.NearestNeighborInterpolator3D
import scalismo.ui.api.TransformationGlyph
import scalismo.mesh.TriangleMesh
import breeze.linalg.DenseVector
import scalismo.statisticalmodel.PointDistributionModel
import scalismo.statisticalmodel.MultivariateNormalDistribution
import breeze.linalg.DenseMatrix
import scalismo.numerics.UniformMeshSampler3D

object ICP {

  /**
    * build a GP Model on a given reference 
    */
  def buildGPModel(
      reference: TriangleMesh[_3D],
      sigma : Double, 
      scaleFactor : Double
  ): LowRankGaussianProcess[_3D, EuclideanVector[_3D]] = {
    val meanFun = (p: Point[_3D]) => EuclideanVector3D(0, 0, 0)
    val mu = Field3D(EuclideanSpace3D, meanFun)

    val cov =
      DiagonalKernel3D(GaussianKernel3D(sigma, scaleFactor), 3)
    val gp = GaussianProcess(mu, cov)
    LowRankGaussianProcess.approximateGPCholesky(
      reference,
      gp,
      relativeTolerance = 1e-1,
      TriangleMeshInterpolator3D()
    )
  }

  /**
   * Finds correspondences between the reference and target mesh, by means
   * of a warped version of the reference mesh (the moving mesh). The 
   * correspondences are only sought on the given point ids
   */
  def attributeCorrespondences(
      movingMesh: TriangleMesh[_3D],
      targetMesh: TriangleMesh[_3D],
      ptIds: Seq[PointId]
  ): Seq[(PointId, Point[_3D])] = {
    ptIds.map { (id: PointId) =>
      val pt = movingMesh.pointSet.point(id)
      val closestPointOnMesh2 = targetMesh.pointSet.findClosestPoint(pt).point
      (id, closestPointOnMesh2)
    }
  }

  /**
   * Condition the model on the given observations (correspondences)
   */
  def fitModel(
      model: PointDistributionModel[_3D, TriangleMesh],
      correspondences: Seq[(PointId, Point[_3D])],
      noiseVariance: Double
  ): TriangleMesh[_3D] = {

    val noiseDist =
      MultivariateNormalDistribution(
        DenseVector.zeros[Double](3),
        DenseMatrix.eye[Double](3) * noiseVariance
      )

    val regressionData =
      correspondences.map(correspondence =>
        (correspondence._1, correspondence._2, noiseDist)
      )
    val posterior = model.posterior(regressionData.toIndexedSeq)
    posterior.mean
  }

  /**
   * Perform the actual ICP iteration
   */
  def nonrigidICP(
      model: PointDistributionModel[_3D, TriangleMesh],
      movingMesh: TriangleMesh[_3D],
      targetMesh: TriangleMesh[_3D],
      ptIds: Seq[PointId],
      numberOfIterations: Int,
      noise: Double
  ): TriangleMesh[_3D] = {
    if (numberOfIterations == 0) movingMesh
    else {
      val correspondences =
        attributeCorrespondences(movingMesh, targetMesh, ptIds)
      val bestFit = fitModel(model, correspondences, noise)

      nonrigidICP(
        model,
        bestFit,
        targetMesh,
        ptIds,
        numberOfIterations - 1,
        noise
      )
    }
  }

  def main(args: Array[String]): Unit = {
    scalismo.initialize()
    val ui = ScalismoUI()

    val faceMesh = MeshIO.readMesh(java.io.File("datasets/lowResPaola.stl")).get
    val gp = buildGPModel(faceMesh, sigma = 100, scaleFactor = 10)

    // we work with a PDM (a GP discretized on the reference mesh) 
    // and not with the raw GP. This makes the intermediate steps easier
    // to follow and visualize
    val model = PointDistributionModel(faceMesh, gp)


    val modelGroup = ui.createGroup("model")
    ui.show(modelGroup, model, "model")

    // to illustrate the principle, we generate a target shape
    // with know coefficients from the model
    // This is called fake data simulation
    val coeffs = DenseVector.zeros[Double](gp.rank)
    coeffs(0) = 2
    coeffs(1) = -2
    val targetMesh = model.instance(coeffs)

    val targetGroup = ui.createGroup("targetGroup")
    ui.show(targetGroup, targetMesh, "targetMesh")

    val ptIds = model.reference.pointSet.pointIds.toIndexedSeq

    val correspondences =
      attributeCorrespondences(model.mean, targetMesh, ptIds)

    val fit =
      nonrigidICP(model.truncate(10), model.mean, targetMesh, ptIds, 100, 1.0)
    val resultGroup = ui.createGroup("results")

    val fitResultView = ui.show(resultGroup, fit, "fit")

  }

}
