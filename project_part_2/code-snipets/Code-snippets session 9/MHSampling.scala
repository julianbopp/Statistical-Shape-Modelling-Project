//> using scala "3.2"
//> using repository "sonatype:snapshots"
//> using dep "ch.unibas.cs.gravis::scalismo-ui:0.91.2"
//> using dep "ch.unibas.cs.gravis::scalismo-plot:0.3-SNAPSHOT"

import breeze.linalg.DenseVector
import breeze.stats.distributions.Rand.FixedSeed.randBasis
import breeze.stats.distributions.Gaussian
import breeze.stats.distributions.MultivariateGaussian
import breeze.linalg.DenseMatrix
import breeze.numerics.step
import scalismo.plot.data.DataFrame
import scalismo.plot.data.DataFrame.Column
import scalismo.plot.plottarget.PlotTargets.plotTargetBrowser

object MHSampling extends App {

    val rng = new scala.util.Random()

    //
    // type definitions to improve clarity of the code
    //
    type Sample = DenseVector[Double]  // The state is represented as a vector
    type Proposal = Sample => Sample    // The proposal function produces a new state from a given state
    type DistributionEvaluator = Sample => Double // The distribution evaluator evaluates a probability of each state

    // 
    // Implementation of the sampler
    //
    def metropolisSampler(p : DistributionEvaluator, q : Proposal)(initialSample : Sample) : Iterator[Sample] = {

        // Simulates one step
        def nextStep(currentSample : Sample) : Sample = {
            
            // propose a new state from the given state
            val proposedSample = q(currentSample) 
            
            // accept based on the ratio of probabilities between 
            // the new and the old state
            val r = rng.nextDouble()
            val alpha  = scala.math.min(1.0, p(proposedSample) / p(currentSample));
            val nextSample = if (r < alpha) proposedSample else currentSample
            nextSample
        }
    
        // create an iterator starting from the initial state
        Iterator.iterate(initialSample)(nextStep)
    }

    //
    // Test program 
    // 

    // target distribution
    val p = (sample : Sample) => {
        val mvn = new MultivariateGaussian(
            DenseVector(3.0, 1.0), 
            DenseMatrix((5.0, 2.0), (2.0, 5.0)))            
        if sample(0) > 2 then 0 else mvn.pdf(sample)
    }
    val stepSize = 0.5
    val q : Proposal = (x : Sample) => new MultivariateGaussian(x, DenseMatrix.eye[Double](2) * stepSize).draw()
    val initialSample = DenseVector(0.0, 0.0)
    
    val samples = metropolisSampler(p, q)(initialSample)
        .drop(1000)  // drop the first n samples to allow for burn-in
        .take(10000)
        .toIndexedSeq
    
    // compute mean and covariance of the samples
    val mean = samples.reduce((sum, s) => sum + s) / samples.size.toDouble
    val covariance = samples.map(s => (s - mean) * (s - mean).t).reduce((sum, s) => sum + s) / samples.size.toDouble
    println(s"mean: $mean")
    println(s"covariance: $covariance")

    // visualize the samples
    DataFrame.fromColumns(Seq(
        Column.ofContinuous(samples.map(s => s(0)), "x"),
        Column.ofContinuous(samples.map(s => s(1)), "y")
    )).plot.scatterPlot("x", "y", "samples").show() 
}