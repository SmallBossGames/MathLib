import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction
import smallBossMathLib.shared.newtonRaphson
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator
import smallBossMathLib.explicitDifferentialEquations.StabilityControlEulerIntegrator
import smallBossMathLib.implicitDifferentialEquations.solveEuler
import kotlin.math.sin

class MyFunction:UnivariateDifferentiableFunction{
    override fun value(t: DerivativeStructure?): DerivativeStructure {
        requireNotNull(t)
        return t.pow(3).multiply(27.0).subtract(t.multiply(3.0)).add(1.0)
    }

    override fun value(x: Double): Double = 27*x*x*x - 3*x + 1
}

class CircleODE(val c: DoubleArray, val omega: Double) : FirstOrderDifferentialEquations {
    override fun computeDerivatives(t: Double, y: DoubleArray?, yDot: DoubleArray?) {
        requireNotNull(y)
        requireNotNull(yDot)

        yDot[0] = omega * (c[1] - y[1])
        yDot[1] = omega * (c[0] - y[0])
    }

    override fun getDimension(): Int = 2
}

fun main() {
    /*val result = newtonRaphson(-0.8, -1.0, 1.0, 8, 1000) {
        x -> Pair(27*x*x*x - 3*x + 1, 27*3*x*x - 3)
    }

    val solver = NewtonRaphsonSolver()
    val func = MyFunction()
    val resultToo = solver.solve(1000, func, -2.0, 3.0)

    println("My lib result: $result, fucking apache: $resultToo");*/

    /*val resultData = implicitEuler(1.0, 3.0, 0.0, 1000){
            x, y -> y.multiply(x).multiply(y)
    }

    println(resultData)*/

    val dp853 = DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10)
    val ode = CircleODE(doubleArrayOf(1.0, 1.0), 0.1)
    val y = doubleArrayOf(0.0, 1.0)

    dp853.integrate(ode, 0.0, y, 16.0, y)

    println("${y[0]}, ${y[1]}")

    /*val result =  solveEuler(0.0, kotlin.math.PI*2, 1000, 0.0){
        t, y -> y.createConstant(t).sin()
    }*/

    val omega = 0.1
    val c = doubleArrayOf(1.0, 1.0)
    val solver = StabilityControlEulerIntegrator(9000, 100)
    val output = DoubleArray(2)

    solver.integrate(0.0, doubleArrayOf(0.0, 1.0), 16.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = omega * (c[1] - inY[1])
        outY[1] = omega * (c[0] - inY[0])
    }

    println(output[0])
    println(output[1])
}