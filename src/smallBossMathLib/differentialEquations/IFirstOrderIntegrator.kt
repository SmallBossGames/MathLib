package smallBossMathLib.differentialEquations

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.FirstOrderIntegrator

interface IFirstOrderIntegrator {
    val maxEvaluations : Int
    val evaluations : Int

    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (t: Double, inY: DoubleArray, outY: DoubleArray) -> Unit
    )
}