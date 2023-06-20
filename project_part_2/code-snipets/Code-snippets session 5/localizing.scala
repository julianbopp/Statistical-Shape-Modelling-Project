//> using scala "3.2.1"
//> using repository "sonatype:snapshots"
//> using lib "ch.unibas.cs.gravis::scalismo-ui:0.92-RC1"

import scalismo.kernels.*
import scalismo.statisticalmodel.{LowRankGaussianProcess, GaussianProcess}
import scalismo.common.*
import scalismo.geometry.*
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.io.{MeshIO, LandmarkIO}

import scalismo.ui.api.ScalismoUI

import scalismo.utils.Random.FixedSeed.randBasis
import scalismo.io.StatisticalModelIO
import scalismo.statisticalmodel.PointDistributionModel


object Localizing extends App {

    scalismo.initialize()
    val ui = ScalismoUI()

    val pdm = StatisticalModelIO.readStatisticalTriangleMeshModel3D(java.io.File("datasets/lowresModel.h5")).get
    val modelGroup1 = ui.createGroup("original model")
    ui.show(modelGroup1, pdm, "original model")


    val pdmGp = pdm.gp.interpolate(TriangleMeshInterpolator3D())
        
    
    
    val gaussianKernel = DiagonalKernel3D(GaussianKernel3D(100, 1), 3)
    
    val localizedGp = GaussianProcess(pdmGp.mean, gaussianKernel * pdmGp.cov)

    val localizedLRGP = LowRankGaussianProcess.approximateGPCholesky(pdm.reference, localizedGp, 1e-2, TriangleMeshInterpolator3D())
    val localizedPDM = PointDistributionModel(pdm.reference, localizedLRGP)

    val modelGroup2 = ui.createGroup("localized model")
    ui.show(modelGroup2, localizedPDM, "localized model")


}   