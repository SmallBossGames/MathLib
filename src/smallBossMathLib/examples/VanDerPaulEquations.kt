package smallBossMathLib.examples

import smallBossMathLib.explicitDifferentialEquations.*
import smallBossMathLib.implicitDifferentialEquations.MK22Integrator
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import java.io.File

fun defaultVanDerPaul(mu: Double, inY: DoubleArray, outF:DoubleArray){
    outF[0] = inY[1]
    outF[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
}

fun pseudoVanDerPaul(mu: Double, inY: DoubleArray, outF:DoubleArray){
    outF[0] = inY[1]
    outF[1] = ((1.0 - inY[0] * inY[1])*inY[1] - inY[0]) / mu
}

fun mk22VdPExample(mu: Double)
{
    val solver = MK22Integrator(
        1e-3,
        0.001,
        20,
        2.0)

    File("VanDerPaul(mu = ${mu}).csv ").bufferedWriter().use { out ->
        out.appendln("t;y1;y2;stiffnessCoefficient")
        solver.addStepHandler { t, y, state, _ ->
            val stiffness = state.jacobiMatrix.evalStiffness()
            val str = "${t};${y[0]};${y[1]};$stiffness"
                .replace('.', ',')
            out.appendln(str);
        }

        val output = doubleArrayOf(2.0, 0.0)
        val rVector = DoubleArray(output.size) { 1e-7 }

        try {
            val result = solver.integrate(0.0, output,20.0, rVector, output)
            {inY, outY ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
            println(result.toString())
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}

fun rk2stVdPExample(mu: Double){
    val solverRK2 = RK23StabilityControlIntegrator(100,0.01)

    File("VanDerPaul(mu = $mu).rk23st.csv ").bufferedWriter().use { out ->
        solverRK2.addStepHandler { t, y, _ ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str)
        }
        //solverRK2.enableStepCountLimit(2000)
        val output = doubleArrayOf(-2.0, 0.0)

        try {
            val result = solverRK2.integrate(0.0, output,20.0, output)
            { _, inY, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
            println("rk2stVdPExample: $result")
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}

fun rk2VdPExample(mu: Double){
    val solverRK2 = RK23Integrator(100,0.01)

    File("VanDerPaul(mu = $mu).rk23.csv ").bufferedWriter().use { out ->
        solverRK2.addStepHandler {t, y, _ ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str)
        }
        //solverRK2.enableStepCountLimit(2000)
        val output = doubleArrayOf(-2.0, 0.0)
        val rVector = DoubleArray(output.size) {1e-7}

        try {
            val result = solverRK2.integrate(0.0, output,20.0, rVector, output)
            { inY: DoubleArray, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
            println("rk2VdPExample: $result")
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}


fun rk4stVdPExample(mu: Double){
    val solverRK2 = RK4StabilityControlIntegrator(100,0.01)

    File("VanDerPaul(mu = $mu).rk45st.csv ").bufferedWriter().use { out ->
        solverRK2.addStepHandler { t, y, _ ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str)
        }
        //solverRK2.enableStepCountLimit(2000)
        val output = doubleArrayOf(-2.0, 0.0)

        try {
            val result = solverRK2.integrate(0.0, output,20.0, output)
            { inY: DoubleArray, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
            println("RK4ST: $result")
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}

fun rkm4VdPExample(mu: Double){
    val solverRK2 = RKM4Integrator(100,0.01)

    File("VanDerPaul(mu = $mu).rkm4.csv ").bufferedWriter().use { out ->
        solverRK2.addStepHandler { t, y, _ ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str)
        }
        //solverRK2.enableStepCountLimit(2000)
        val output = doubleArrayOf(-2.0, 0.0)
        val rVector = DoubleArray(output.size) {1e-7}

        try {
            val result = solverRK2.integrate(0.0, output,20.0, rVector, output)
            { inY: DoubleArray, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
            println("RKM4: $result")

        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}


fun mk22VdPAlternateExample(p: Double)
{
    val solver = MK22Integrator(
        1e-3,
        0.001,
        0,
        0.0)

    File("VanDerPaul(p = ${p}).mk22.csv ").bufferedWriter().use { out ->
        out.appendln("t;y1;y2;stiffnessCoefficient")
        solver.addStepHandler { t, y, state, _ ->
            val stiffness = state.jacobiMatrix.evalStiffness()
            val str = "${t};${y[0]};${y[1]};$stiffness"
                .replace('.', ',')

            out.appendln(str);
        }
        solver.enableStepCountLimit(20000)

        val output = doubleArrayOf(2.0, 0.0)
        val rVector = DoubleArray(output.size) { 1e-7 }

        try {
            solver.integrate(0.0, output,1.0, rVector, output)
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

fun eulerVdPAlternateExample(p: Double)
{
    val solver = EulerIntegrator(10000, 0.01)

    File("VanDerPaul(p = ${p}).csv ").bufferedWriter().use { out ->
        out.appendln("t;y1;y2")
        solver.addStepHandler { t, y, _ ->
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

fun eulerVdPExample(mu: Double)
{
    val solver = EulerIntegrator(
        10000,
        0.01)

    File("VanDerPaul(mu = ${mu}).euler.csv ").bufferedWriter().use { out ->
        out.appendln("t;y1;y2")
        solver.addStepHandler { t, y, _ ->
            val str = "${t};${y[0]};${y[1]}".replace('.', ',')
            out.appendln(str);
        }
        solver.enableStepCountLimit(20000)

        val output = doubleArrayOf(2.0, 0.0)

        try {
            val result = solver.integrate(0.0, output,20.0, output)
            {inY: DoubleArray, outY: DoubleArray ->
                outY[0] = inY[1]
                outY[1] = mu * (1 - inY[0] * inY[0]) * inY[1] - inY[0]
            }
            println("eulerVdPExample: $result")
        } catch (ex: ExceedingLimitEvaluationsException){
            println(ex.message)
        } catch (ex: ExceedingLimitStepsException){
            println(ex.message)
        }
    }
}