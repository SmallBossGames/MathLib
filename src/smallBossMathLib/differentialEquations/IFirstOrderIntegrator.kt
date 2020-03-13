package smallBossMathLib.differentialEquations

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.FirstOrderIntegrator

interface IFirstOrderIntegrator {
    val maxEvaluationCount : Int
    val evaluations : Int
    val accuracy : Double

    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (t: Double, inY: DoubleArray, outY: DoubleArray) -> Unit
    )

    fun addStepHandler(handler: (t: Double, y: DoubleArray) -> Unit)
}