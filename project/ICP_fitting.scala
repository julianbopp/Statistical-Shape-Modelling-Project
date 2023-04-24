//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using lib "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import scalismo.geometry._
import scalismo.common._
import scalismo.mesh._
import scalismo.statisticalmodel.MultivariateNormalDistribution
import scalismo.numerics.UniformMeshSampler3D
import scalismo.io.{MeshIO, StatisticalModelIO, LandmarkIO}

import scalismo.ui.api._

import breeze.linalg.{DenseMatrix, DenseVector}

object ICP_fitting extends App {
    
    scalismo.initialize()
    implicit val rng: scalismo.utils.Random = scalismo.utils.Random(42)

    val ui = ScalismoUI()

    val model = StatisticalModelIO.readStatisticalTriangleMeshModel3D(new java.io.File("gpmodel.h5")).get
    
    val targetGroup = ui.createGroup("targetGroup")

    val modelGroup = ui.createGroup("modelGroup")
    val modelView = ui.show(modelGroup, model, "model")
    val resultGroup = ui.createGroup("results")

    val sampler = UniformMeshSampler3D(model.reference, numberOfPoints = 5000)
    val points: Seq[Point[_3D]] = sampler.sample().map(pointWithProbability => pointWithProbability._1) // we only want the points
    val ptIds = points.map(point => model.reference.pointSet.findClosestPoint(point).id)

    val referenceLandmarks = LandmarkIO
      .readLandmarksJson3D(
        new java.io.File("project-data/reference-landmarks/reference.json")
      )
      .get

    val refLandmarkPoints = referenceLandmarks.map(_.point)
    val refLandmarkPointIds = refLandmarkPoints.map(point => model.reference.pointSet.findClosestPoint(point))
    print(refLandmarkPointIds)

    val numOfTargets = 47

    for (i <- 0 until numOfTargets) {
        val targetMesh = MeshIO.readMesh(new java.io.File(s"results/aligned/$i.stl")).get
        ui.show(targetGroup, targetMesh, s"targetMesh$i")

        // Find for each point of interest the closest point on the target
        def attributeCorrespondences(movingMesh: TriangleMesh[_3D], ptIds: Seq[PointId]): Seq[(PointId, Point[_3D])] = {
            ptIds.map { (id: PointId) => val pt = movingMesh.pointSet.point(id)
                val closestPointOnMesh2 = targetMesh.pointSet.findClosestPoint(pt).point
                (id, closestPointOnMesh2)
            }
        }

        val correspondences = attributeCorrespondences(model.mean, ptIds)

        val littleNoise = MultivariateNormalDistribution(DenseVector.zeros[Double](3), DenseMatrix.eye[Double](3))

        def fitModel(correspondences: Seq[(PointId, Point[_3D])]) : TriangleMesh[_3D] = {
            val regressionData = correspondences.map(correspondence =>
            (correspondence._1, correspondence._2, littleNoise)
            )
            val posterior = model.posterior(regressionData.toIndexedSeq)
            posterior.mean
        }

        val fit = fitModel(correspondences)
        ui.show(resultGroup, fit, s"fit $i")

        def nonrigidICP(
            movingMesh: TriangleMesh[_3D],
            ptIds: Seq[PointId],
            numberOfIterations: Int
        ): TriangleMesh[_3D] = {
        if (numberOfIterations == 0) movingMesh
        else {
            val correspondences = attributeCorrespondences(movingMesh, ptIds)
            val transformed = fitModel(correspondences)

            nonrigidICP(transformed, ptIds, numberOfIterations - 1)
        }
        }
        val finalFit = nonrigidICP(model.mean, ptIds, 20)
        MeshIO
            .writeMesh(finalFit, new java.io.File(s"results/fitted/fitted_$i.stl"))
            .get

        ui.show(resultGroup, finalFit, s"final fit $i") 
    }
}
