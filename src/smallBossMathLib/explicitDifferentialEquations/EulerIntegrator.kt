package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.zeroSafetyNorm
import kotlin.math.sqrt


class EulerIntegrator(val defaultStep: Double,
                      val accuracy: Double) : ExplicitIntegrator() {
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (y: DoubleArray, f: DoubleArray) -> Unit
    ) : IExplicitMethodStatistic {
        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fNextBuffer = DoubleArray(y0.size)

        val yNextBuffer = DoubleArray(y0.size)

        val vectorBuffer = DoubleArray(y0.size)

        var time = t0
        var step = defaultStep

        val endTime = time+t

        val statistic = ExplicitMethodStatistic(
            stepsCount = 0,
            evaluationsCount = 0,
            returnsCount = 0
        )

        val state = ExplicitMethodStepState(
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step)
        )

        executeStepHandlers(time, outY, state, statistic)

        equations(outY, fCurrentBuffer)

        while (time < endTime){
            step = normalizeStep(step, time, endTime)

            for (i in yNextBuffer.indices){
                yNextBuffer[i] = outY[i] + step*fCurrentBuffer[i]
            }

            equations(yNextBuffer, fNextBuffer)

            statistic.evaluationsCount++

            for (i in vectorBuffer.indices){
                vectorBuffer[i] = fNextBuffer[i] - fCurrentBuffer[i]
            }

            val errNorm = 0.5 * step * zeroSafetyNorm(vectorBuffer, outY, 1e-7)
            val q = sqrt(accuracy / errNorm)

            if(q < 1.0 && !state.isLowStepSizeReached && !state.isHighStepSizeReached){
                step = q * step / 1.1

                state.isLowStepSizeReached = isLowStepSizeReached(step)
                state.isHighStepSizeReached = isHighStepSizeReached(step)
                statistic.returnsCount++
            }
            else {
                for (i in outY.indices){
                    outY[i] = yNextBuffer[i]
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