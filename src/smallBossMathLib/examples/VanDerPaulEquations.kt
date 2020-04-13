package smallBossMathLib.examples

import smallBossMathLib.explicitDifferentialEquations.RK23StabilityControlIntegrator
import smallBossMathLib.implicitDifferentialEquations.MK22Integrator
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import java.io.File

fun rungeKuttaSecondOrderExample(p: Double)
{
    val builder = StringBuilder()
    val solverRK2 = RK23StabilityControlIntegrator(100,0.01)
    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}
    solverRK2.enableStepCountLimit(10000)

    val output = doubleArrayOf(2.0, 0.0)

    solverRK2.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = ((1.0 - inY[0] * inY[1])*inY[1] - inY[0]) / p
    }

    val writingText = builder.replace(Regex("[.]"), ",")

    File("VanDerPaul(p = ${p}).csv ").writeText(writingText)
}

fun mk22VdPExample(mu: Double)
{
    val solver = MK22Integrator(
        10000,
        0,
        0.0,
        0.001,
        0.0,
        Double.POSITIVE_INFINITY)

    File("VanDerPaul(mu = ${mu}).csv ").bufferedWriter().use { out ->
        out.appendln("t;y1;y2")
        solver.addStepHandler { t, y ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str);
        }
        solver.enableStepCountLimit(20000)

        val output = doubleArrayOf(2.0, 0.0)

        try {
            solver.integrate(0.0, output,20.0, output)
            {inY: DoubleArray, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}

fun rk2VdPExample(mu: Double){
    val builder = StringBuilder()
    val solverRK2 = RK23StabilityControlIntegrator(100,0.001)
    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}
    solverRK2.enableStepCountLimit(2000)

    val output = doubleArrayOf(-2.0, 0.0)

    solverRK2.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
    }

    for (i in output)
        println(i)
    //val writingText = builder.replace(Regex("[.]"), ",")

    //File("data(${mu}).csv ").writeText(writingText)
}

fun mk22VdPAlternateExample(p: Double)
{
    val solver = MK22Integrator(
        10000,
        0,
        0.0,
        0.01,
        0.0,
        Double.POSITIVE_INFINITY)

    File("VanDerPaul(p = ${p}).csv ").bufferedWriter().use { out ->
        out.appendln("t;y1;y2")
        solver.addStepHandler { t, y ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str);
        }
        solver.enableStepCountLimit(20000)

        val output = doubleArrayOf(2.0, 0.0)

        try {
            solver.integrate(0.0, output,10.0, output)
            {inY: DoubleArray, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = ((1.0 - inY[0] * inY[1])*inY[1] - inY[0]) / p
            }
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}