package smallBossMathLib.examples

import smallBossMathLib.explicitDifferentialEquations.StabilityControlSecondOrderIntegrator
import smallBossMathLib.implicitDifferentialEquations.MK22Integrator
import java.io.File

fun rungeKuttaSecondOrderExample(p: Double)
{
    val builder = StringBuilder()
    val solverRK2 = StabilityControlSecondOrderIntegrator(100,0.01)
    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}
    solverRK2.enableEvaluationCountCheck(10000)

    val output = doubleArrayOf(2.0, 0.0)

    solverRK2.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = ((1.0 - inY[0] * inY[1])*inY[1] - inY[0]) / p
    }

    val writingText = builder.replace(Regex("[.]"), ",")

    File("VanDerPaul(p = ${p}).csv ").writeText(writingText)
}

fun mk22VdPExample(mu: Double){
    val builder = StringBuilder()
    val solverRK2 = MK22Integrator(
        100,
        0,
        0.0,
        0.001,
        0.0,
        Double.POSITIVE_INFINITY)

    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}
    //solverRK2.enableEvaluationCountCheck(2000)

    val output = doubleArrayOf(-2.0, 0.0)

    solverRK2.integrate(0.0, output,20.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
    }

    //for (i in output)
    //    println(i)
    val writingText = builder.replace(Regex("[.]"), ",")

    File("mk22test(${mu}).csv ").writeText(writingText)
}

fun rk2VdPExample(mu: Double){
    val builder = StringBuilder()
    val solverRK2 = StabilityControlSecondOrderIntegrator(100,0.001)
    solverRK2.addStepHandler { t, y ->  builder.append("${y[0]};${y[1]} \n");}
    solverRK2.enableEvaluationCountCheck(2000)

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
    val builder = StringBuilder()
    val solver = MK22Integrator(
        10000,
        0,
        0.0,
        0.001,
        0.0,
        Double.POSITIVE_INFINITY)

    solver.addStepHandler { t, y ->  builder.appendln("${t};${y[0]};${y[1]}");}
    solver.enableEvaluationCountCheck(20000)
    builder.appendln("t;y1;y2")

    val output = doubleArrayOf(2.0, 0.0)

    solver.integrate(0.0, output,10.0, output)
    { t: Double, inY: DoubleArray, outY: DoubleArray ->
        outY[0] = inY[1]
        outY[1] = ((1.0 - inY[0] * inY[1])*inY[1] - inY[0]) / p
    }

    val writingText = builder.replace(Regex("[.]"), ",")

    File("VanDerPaul(p = ${p}).csv ").writeText(writingText)
}