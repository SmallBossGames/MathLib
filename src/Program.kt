import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import smallBossMathLib.explicitDifferentialEquations.StabilityControlFourthOrderIntegrator
import smallBossMathLib.explicitDifferentialEquations.StabilityControlSecondOrderIntegrator
import smallBossMathLib.examples.*
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
    //Let's try with Van-der-Paul

    /*testRungeKuttaWithStepHandlerSecondOrder(8.0)
    testRungeKuttaWithStepHandlerSecondOrder(6.0)
    testRungeKuttaWithStepHandlerSecondOrder(4.0)
    testRungeKuttaWithStepHandlerSecondOrder(3.0)
    testRungeKuttaWithStepHandlerSecondOrder(2.5)
    testRungeKuttaWithStepHandlerSecondOrder(2.0)
    testRungeKuttaWithStepHandlerSecondOrder(1.0)
    testRungeKuttaWithStepHandlerSecondOrder(0.01)*/

    //testRungeKuttaWithStepHandlerFourthOrder(1000.0)

    //ringModualtorExample()

    //rungeKuttaSecondOrderExample(3e-2)

    luTest()


    //rk2VdPExample(6.0)


    ringModulatorMK22Example()
    mk22VdPExample(6.0)
    //mk22Other2Test()
    //mk22Other3Test()

    mk22VdPAlternateExample(3.0e-6)
    mk22VdPAlternateExample(3.0e-5)
    mk22VdPAlternateExample(3.0e-4)
    mk22VdPAlternateExample(3.0e-3)
    mk22VdPAlternateExample(3.0e-2)
}

fun testRungeKuttaSecondOrder(mu: Double)
{
    val solverRK2 = StabilityControlSecondOrderIntegrator(100, 0.01)

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
    val solverRK2 = StabilityControlSecondOrderIntegrator(100,0.01)
    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}
    //solverRK2.enableEvaluationCountCheck(2000)

    val output = doubleArrayOf(-2.0, 0.0)

    solverRK2.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
    }

    val writingText = builder.replace(Regex("[.]"), ",")

    File("data(${mu}).csv ").writeText(writingText)
}

fun testRungeKuttaWithStepHandlerFourthOrder(mu: Double)
{
    val builder = StringBuilder()
    val solverRK4 = StabilityControlFourthOrderIntegrator( 100,0.01)
    solverRK4.addStepHandler { _, y ->  builder.append("${y[0]};${y[1]} \n");}

    val output = doubleArrayOf(-2.0, 0.0)

    solverRK4.integrate(0.0, output,2000.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
    }

    val writingText = builder.replace(Regex("[.]"), ",")

    File("data4(${mu}).csv ").writeText(writingText)
}