import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator
import org.apache.commons.math3.ode.nonstiff.MidpointIntegrator
import smallBossMathLib.explicitDifferentialEquations.StabilityControlFourthOrderIntegrator
import smallBossMathLib.explicitDifferentialEquations.StabilityControlSecondOrderIntegrator
import java.io.File

class MyFunction:UnivariateDifferentiableFunction{
    override fun value(t: DerivativeStructure?): DerivativeStructure {
        requireNotNull(t)
        return t.pow(3).multiply(27.0).subtract(t.multiply(3.0)).add(1.0)
    }

    override fun value(x: Double): Double = 27*x*x*x - 3*x + 1
}

class VanDerPol(private val mu: Double) : FirstOrderDifferentialEquations {
    override fun computeDerivatives(t: Double, y: DoubleArray?, yDot: DoubleArray?) {
        requireNotNull(y)
        requireNotNull(yDot)

        yDot[0] = y[1]
        yDot[1] = mu * (1 - y[0] * y[0]) * y[1] - y[0]
    }

    override fun getDimension(): Int = 2
}

class CircleODE(private val c: DoubleArray, private val omega: Double) : FirstOrderDifferentialEquations {
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

    //val dp853 = DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10)

    /*val ode = CircleODE(doubleArrayOf(1.0, 1.0), 0.1)
    val y = doubleArrayOf(0.0, 1.0)

    dp853.integrate(ode, 0.0, y, 16.0, y)

    println("${y[0]}, ${y[1]}")*/


   // val y = doubleArrayOf(2.0, 0.0)

    //dp853.integrate(vanDerPol, 0.0, y, 3000.0, y)
    //println("${y[0]}, ${y[1]}")



    /*val omega = 0.1
    val c = doubleArrayOf(1.0, 1.0)


    solver.integrate(0.0, doubleArrayOf(0.0, 1.0), 16.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = omega * (c[1] - inY[1])
        outY[1] = omega * (c[0] - inY[0])
    }

    println(output[0])
    println(output[1])*/

    //Let's try with Van-der-Paul

    /*testRungeKuttaWithStepHandlerSecondOrder(8.0)
    testRungeKuttaWithStepHandlerSecondOrder(6.0)
    testRungeKuttaWithStepHandlerSecondOrder(4.0)
    testRungeKuttaWithStepHandlerSecondOrder(3.0)
    testRungeKuttaWithStepHandlerSecondOrder(2.5)
    testRungeKuttaWithStepHandlerSecondOrder(2.0)
    testRungeKuttaWithStepHandlerSecondOrder(1.0)
    testRungeKuttaWithStepHandlerSecondOrder(0.01)*/

    testRungeKuttaWithStepHandlerFourthOrder(8.0)


}

fun testRungeKuttaSecondOrder(mu: Double)
{
    val solverRK2 = StabilityControlSecondOrderIntegrator(9000, 100, 0.01)

    for(i in 0..400)
    {
        val output = doubleArrayOf(-2.0, 0.0)

        solverRK2.integrate(0.0, output,0.05 * i, output)
        { t: Double, inY: DoubleArray, outY: DoubleArray ->
            outY[0] = inY[1]
            outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
        }

        print("(${output[0]};${output[1]}) ")
    }
}

fun testRungeKuttaWithStepHandlerSecondOrder(mu: Double)
{
    val builder = StringBuilder()
    val solverRK2 = StabilityControlSecondOrderIntegrator(9000, 100,0.01)
    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}

    val output = doubleArrayOf(-2.0, 0.0)

    solverRK2.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
    }

    File("data(${mu}).csv ").writeText(builder.toString())
}

fun testRungeKuttaWithStepHandlerFourthOrder(mu: Double)
{
    val builder = StringBuilder()
    val solverRK4 = StabilityControlFourthOrderIntegrator(9000, 100,0.01)
    solverRK4.addStepHandler { _, y ->  builder.append("${y[0]};${y[1]} \n");}

    val output = doubleArrayOf(-2.0, 0.0)

    solverRK4.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
    }

    File("data4(${mu}).csv ").writeText(builder.toString())
}