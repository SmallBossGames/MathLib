import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import smallBossMathLib.examples.*
import kotlin.system.measureTimeMillis

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
    //ringModualtorExample()

    //rungeKuttaSecondOrderExample(3e-2)

    //luTest()

    //kroneckerMultiplyTest1()

    //kroneckerMultiplyTest2()

    //rk2VdPExample(6.0)
    
    //ringModulatorMK22Example()
    //ringModulatorRK23ST()
    //ringModulatorImplicitEulerExample()

    //mk22Other2Test()
    //mk22Other3Test()

    //mk22VdPAlternateExample(1e-6)
    //mk22VdPAlternateExample(3.0e-5)
    //mk22VdPAlternateExample(3.0e-4)
    //mk22VdPAlternateExample(3.0e-3)
    //mk22VdPAlternateExample(0.75e-2)
    //mk22VdPAlternateExample(1.5e-2)


    //mk22VdPAlternateExample(3e-2)

    //mk22VdPExample(0.01)
    //mk22VdPExample(6.0)

    //eulerVdPExample(6.0)

    //eulerVdPAlternateExample(1e-2)

    //implicitEulerVdPExample(6.0)

    //rk2stVdPExample(6.0)

    //rk2VdPExample(6.0)

    //rk4stVdPExample(6.0)

    //rkm4VdPExample(6.0)

    //radau5Order3VdPExample(6.0)

    //ringModulatorRadau5Example()

    //impEulerAbsoluteAccuracyTest()
    //mk22VdPAbsoluteAccuracyTest()
    //eulerAbsoluteAccuracyTest()
    //rk23AbsoluteAccuracyTest()
    //rkm45AbsoluteAccuracyTest()
    // rk23stAbsoluteAccuracyTest()
    rkm45AbsoluteAccuracyTest()


    //eulerBreakpointTest()
}