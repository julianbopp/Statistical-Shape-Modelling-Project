//> using scala "3.1.2"
//> using lib "ch.unibas.cs.gravis::scalismo-ui:0.91.0"
import scalismo.geometry._
import scalismo.common._
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.mesh._
import scalismo.io.{StatisticalModelIO, MeshIO}
import scalismo.statisticalmodel._
import scalismo.numerics.UniformMeshSampler3D
import scalismo.kernels._

import scalismo.ui.api._

import breeze.linalg.{DenseMatrix, DenseVector}

object GPmodel extends App {

    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()


    // Load reference mesh and display in separate group
    val referenceMesh = MeshIO.readMesh(new java.io.File("project-data/reference-mesh/reference.stl")).get
    val modelGroup = ui.createGroup("gp-model")
    val referenceView = ui.show(modelGroup, referenceMesh, "reference")

    // Define mean deformation as zero mean (i.e. we assume reference shape to be mean shape)
    val zeroMean = Field(EuclideanSpace3D, (pt:Point[_3D]) => EuclideanVector3D(0,0,0))

    // Define Gaussian Kernel
    val scalarValuedGaussianKernel : PDKernel[_3D]= GaussianKernel3D(sigma = 100.0)
    val matrixValuedGaussianKernel = DiagonalKernel3D(scalarValuedGaussianKernel, 3) 
}
