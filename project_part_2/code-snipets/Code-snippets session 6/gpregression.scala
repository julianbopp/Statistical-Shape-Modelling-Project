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

import scalismo.utils.Random.FixedSeed.randBasis
import scalismo.common.interpolation.NearestNeighborInterpolator3D
import scalismo.ui.api.TransformationGlyph
import breeze.stats.distributions.MultivariateGaussian
import breeze.linalg.DenseMatrix
import breeze.linalg.DenseVector

import breeze.stats.distributions.Rand.FixedSeed.randBasis
import scalismo.statisticalmodel.MultivariateNormalDistribution


object GPRegression extends App {

    scalismo.initialize()
    val ui = ScalismoUI()

    val meanFun = (p : Point[_3D]) => EuclideanVector3D(0, 0, 0)
    val mu = Field3D(EuclideanSpace3D, meanFun)

    val cov = DiagonalKernel3D(GaussianKernel3D(sigma = 100, scaleFactor = 10), 3)
    val gp = GaussianProcess(mu, cov)

    val faceMesh = MeshIO.readMesh(java.io.File("datasets/lowResPaola.stl")).get
    ui.show(faceMesh, "reference mesh")

    val landmarks = LandmarkIO.readLandmarksJson3D(java.io.File("datasets/gp-regression-lms.json")).get
    ui.show(landmarks,  "landmarks")

    val sigma = 1.0
    val sigma2 = sigma * sigma

    val noiseDist = MultivariateNormalDistribution(DenseVector.zeros[Double](3), DenseMatrix.eye[Double](3) * sigma2)
    val data = IndexedSeq((landmarks(0).point, EuclideanVector(0, 0, 0), noiseDist), (landmarks(1).point, EuclideanVector(0, 0, 0), noiseDist))

    val posteriorGP = GaussianProcess.regression(gp, data)
    val discreteGP = posteriorGP.discretize(faceMesh)

    ui.show(discreteGP.sample(), "sample")

}   
