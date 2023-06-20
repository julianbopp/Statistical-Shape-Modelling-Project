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


object GPs extends App {

    scalismo.initialize()
    val ui = ScalismoUI()
    
    val meanFun = (p : Point[_3D]) => EuclideanVector3D(0.0, 0.0, 0.0)
    
    val mu = Field3D(EuclideanSpace3D, meanFun)
    val cov = DiagonalKernel3D(GaussianKernel3D(100, 10), 3)

    val gp = GaussianProcess(mu, cov)


    val faceMesh = MeshIO.readMesh(java.io.File("datasets/lowResPaola.stl")).get
    ui.show(faceMesh, "reference face")


    // discretize then sample
    val discreteGP =  gp.discretize(faceMesh)
    val sampleDfDiscrete = discreteGP.sample()

    ui.show(sampleDfDiscrete, "sample deformation field")

    // warp a mesh. We need to make the deformation field continuous, otherwise we cannot
    // guarantee that it is defined everywhere
    val sampleDfCont = sampleDfDiscrete.interpolate(TriangleMeshInterpolator3D())
    val warpedMesh = faceMesh.transform(p => p + sampleDfCont(p))
    ui.show(warpedMesh, " sampled face mesh")


    // approximate to finite rank then discretize
    val lowrankGP = LowRankGaussianProcess.approximateGPCholesky(faceMesh, gp, 1e-5, TriangleMeshInterpolator3D())
    val sampleDfCont2 = lowrankGP.sample()

    val sampleDfDiscrete2 = sampleDfCont2.discretize(faceMesh, outsideValue = EuclideanVector3D(0, 0, 0))
    ui.show(sampleDfDiscrete2, " sample deformation field 2")

    // alternatively, we discretize the gp  
    val discreteLowrankGP = lowrankGP.discretize(faceMesh)
    val sampleDfDiscrete3 = discreteLowrankGP.sample()


    
    // extract interesting measurements
    val landmarks = LandmarkIO.readLandmarksJson3D(java.io.File("datasets/intereye-lms.json")).get
    ui.show(landmarks, "landmarks")

    val p1 = landmarks(0).point
    val p2 = landmarks(1).point

    val nSamples = 100
    val distances = for i <- 0 until nSamples yield 
        val sampleDF = lowrankGP.sample()
        val p1transformed = p1 + sampleDF(p1)
        val p2transformed = p2 + sampleDF(p2)
        (p2transformed - p1transformed).norm

    // visualizing the measurements
    val dataFrame = DataFrame.fromColumns(Seq(
        Column.ofContinuous(distances, "inter-eye-distances")
    )
    ).plot.histogram("inter-eye-distances", "inter-eye distances").show()

}   