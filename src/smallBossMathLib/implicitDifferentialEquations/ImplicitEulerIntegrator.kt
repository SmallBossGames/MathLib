package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitEvaluationsException
import smallBossMathLib.implicitDifferentialEquations.exceptions.ExceedingLimitStepsException
import smallBossMathLib.shared.Matrix2D
import smallBossMathLib.shared.NewtonRaphsonSolver
import smallBossMathLib.shared.StationaryODE
import smallBossMathLib.shared.zeroSafetyNorm
import kotlin.math.sqrt

class ImplicitEulerIntegrator(
    val defaultStep: Double,
    val accuracy: Double
) : ImplicitIntegrator() {
    @Throws(ExceedingLimitStepsException::class, ExceedingLimitEvaluationsException::class)
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        rVector: DoubleArray,
        outY: DoubleArray,
        equations: StationaryODE
    ) : IImplicitMethodStatistic {
        if (y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        val endTime = t0 + t

        var fCurrentBuffer = DoubleArray(y0.size)
        var fNextBuffer = DoubleArray(y0.size)

        val nextY = DoubleArray(y0.size)
        val vectorBuffer1 = DoubleArray(y0.size)

        var step = defaultStep
        var time = t0

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

        equations(y0, fCurrentBuffer)

        while (time < endTime){
            step = normalizeStep(step, time, endTime)

            newtonSolver.solve(outY, nextY) { y, out ->
                equations(y, fNextBuffer)
                for (i in out.indices){
                    out[i] = outY[i] + step*fNextBuffer[i] - y[i]
                }
            }

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = fNextBuffer[i] - fCurrentBuffer[i]
            }

            val errNorm = 0.5 * step * zeroSafetyNorm(vectorBuffer1, outY, rVector)
            val q = sqrt(accuracy / errNorm)

            if(q < 1.0 && !state.isLowStepSizeReached && !state.isHighStepSizeReached){
                step = q * step / 1.1

                state.isLowStepSizeReached = isLowStepSizeReached(step)
                state.isHighStepSizeReached = isHighStepSizeReached(step)
                statistic.returnsCount++
            }
            else {
                for (i in outY.indices){
                    outY[i] = nextY[i]
                }

                time += step

                statistic.stepsCount++

                executeStepHandlers(time, outY, state, statistic)

                step = q * step / 1.1

                state.isLowStepSizeReached = isLowStepSizeReached(step)
                state.isHighStepSizeReached = isHighStepSizeReached(step)

                val temp = fCurrentBuffer
                fCurrentBuffer = fNextBuffer
                fNextBuffer = temp
            }
        }

        return statistic
    }
}