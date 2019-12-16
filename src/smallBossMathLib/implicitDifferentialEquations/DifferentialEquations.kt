package smallBossMathLib.implicitDifferentialEquations

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver
import smallBossMathLib.shared.newtonRaphson
import kotlin.math.abs;

//В общем случае здесь
//yn1 = yn + h * (b1*k1 + b2*k2 + .. bs*ks)
//ki = f(tn + ci * h, yn + h * sum( from j=1 to s aij*kj))

fun radau2A(from: Double, to: Double, y0: Double = 0.0, step: Double, f: (t:Double, y:Double) -> Double):Double {
    val stepCount = abs((to - from)/step).toInt()
    var y = y0

    for (a in 0 until stepCount) {
        val t = a * step

        //Here we have to found k-values with Newton-Raphson method

        //val func = smallBossMathLib.shared.UnivariateDifferentiableKotlinFunction({})

        //val k1 = f(t + 1/3, y + step * (5/12*))

    }

    TODO("Write tests and realisation")
}

data class EulerInnerParameters(var y: Double, var h: Double, var t: Double);

class EulerNewYCalculationFunc(var params: EulerInnerParameters,
                               val f: (t: Double, y: DerivativeStructure) -> DerivativeStructure)
    : UnivariateDifferentiableFunction {
    override fun value(y1: DerivativeStructure?): DerivativeStructure {
        requireNotNull(y1)
        return (f(params.t, y1).multiply(params.h)).add(params.y).subtract(y1)
    }

    override fun value(x: Double): Double {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

}

fun solveEuler(from: Double, to: Double, stepCount: Int, y0: Double,
               f: (t: Double, y: DerivativeStructure) -> DerivativeStructure): Double {
    val solver = NewtonRaphsonSolver()
    val params = EulerInnerParameters(y0, (to - from) / stepCount, 0.0)
    val mainFunc = EulerNewYCalculationFunc(params, f)

    for (i in 0 until stepCount)
    {
        params.t = i * params.h
        params.y = solver.solve(1000, mainFunc, 0.0, 1000.0)
    }
    return params.y
}