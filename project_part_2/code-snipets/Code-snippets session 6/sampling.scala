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
import breeze.stats.distributions.MultivariateGaussian
import breeze.linalg.DenseMatrix
import breeze.linalg.DenseVector

import breeze.stats.distributions.Rand.FixedSeed.randBasis
import scalismo.utils.Random.implicits.randomGenerator


object Sampling extends App {

    scalismo.initialize()
    val ui = ScalismoUI()

    val meanFun = (p : Point[_3D]) => EuclideanVector3D(0, 0, 0)
    val mu = Field3D(EuclideanSpace3D, meanFun)

    val cov = DiagonalKernel3D(GaussianKernel3D(sigma = 100, scaleFactor = 10), 3)
    val gp = GaussianProcess(mu, cov)

    val faceMesh = MeshIO.readMesh(java.io.File("datasets/lowResPaola.stl")).get
    ui.show(faceMesh, "reference mesh")

    val lowrankGP = LowRankGaussianProcess.approximateGPCholesky(faceMesh, gp, relativeTolerance = 1e-1, TriangleMeshInterpolator3D())


    // sample from prior
    val sampleDeformationPrior = lowrankGP.sample()
    val sampleDeformationPriorDiscrete = sampleDeformationPrior.discretize(faceMesh, outsideValue = EuclideanVector3D.zero)

    ui.show(sampleDeformationPriorDiscrete, "sampled deformation prior")

    // sample from likelihood
    val sigma = 1.0
    val sigma2 = sigma * sigma        

    val noiseDistribution = MultivariateGaussian(DenseVector.zeros[Double](3), DenseMatrix.eye[Double](3) * sigma2)

    val sampleDeformationData = sampleDeformationPrior.andThen(u => u +  EuclideanVector.fromBreezeVector[_3D](noiseDistribution.sample()))

    val sampleDeformationDataDiscrete = sampleDeformationData.discretize(faceMesh, outsideValue = EuclideanVector3D.zero)

    ui.show(sampleDeformationDataDiscrete, "sampled deformation data")
    
    val sampleMesh = faceMesh.transform(p => p + sampleDeformationData(p))

    ui.show(sampleMesh, "sampled mesh")

}   
