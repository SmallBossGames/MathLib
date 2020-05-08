package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import smallBossMathLib.shared.NewtonRaphsonSolver
import smallBossMathLib.shared.StationaryODE

class Radau5Order3Integrator(val accuracy: Double) : ImplicitIntegrator() {
    @Throws(ExceedingLimitStepsException::class, ExceedingLimitEvaluationsException::class)
    fun integrate(
        startTime: Double,
        endTime: Double,
        defaultStepSize: Double,
        y0: DoubleArray,
        rVector: DoubleArray,
        outY: DoubleArray,
        equations: StationaryODE
    ) : IImplicitMethodStatistic {
        if (y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        val vectorBuffer1 = DoubleArray(y0.size)
        val vectorBuffer2 = DoubleArray(y0.size)

        val k = Array(2){DoubleArray(y0.size)}

        var step = defaultStepSize
        var time = startTime

        val newtonSolver = NewtonRaphsonSolver(y0.size, accuracy, rVector)

        val statistic = ImplicitMethodStatistic(
            stepsCount = 0,
            evaluationsCount = 0,
            jacobiEvaluationsCount = 0,
            returnsCount = 0
        )

        val state = ImplicitMethodStepState(
            jacobiMatrix = newtonSolver.jacobiMatrix,
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step),
            freezeJacobiStepsCount = 0
        )

        while (time < endTime){
            step = normalizeStep(step, time, endTime)

            equations(outY, vectorBuffer1)
            statistic.evaluationsCount++
            for (i in vectorBuffer1.indices) {
                k[0][i] = step * vectorBuffer1[i]
                vectorBuffer2[i] = outY[i] + 1.0 / 3.0 * k[0][i]
            }

            equations(vectorBuffer2, vectorBuffer1)
            statistic.evaluationsCount++
            for (i in vectorBuffer1.indices) {
                k[0][i] = step * vectorBuffer1[i]
                vectorBuffer1[i] = outY[i] + 1.0 / 3.0 * k[0][i]
            }
        }

        return statistic
    }
}