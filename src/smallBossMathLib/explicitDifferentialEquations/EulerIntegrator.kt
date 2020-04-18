package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.Integrator
import smallBossMathLib.shared.StepInfo
import smallBossMathLib.shared.zeroSafetyNorm
import kotlin.math.sqrt


class EulerIntegrator(val startEvaluations: Int,
                      val accuracy: Double) : ExplicitIntegrator() {
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (y: DoubleArray, f: DoubleArray) -> Unit
    ) : IExplicitMethodStepInfo {
        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fNextBuffer = DoubleArray(y0.size)

        val yNextBuffer = DoubleArray(y0.size)

        val vectorBuffer = DoubleArray(y0.size)

        var time = t0
        var step = t / startEvaluations

        val endTime = time+t

        val stepInfo = ExplicitMethodStepInfo(
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step),
            stepsCount = 0,
            evaluationsCount = 0,
            returnsCount = 0
        )

        executeStepHandlers(time, outY, stepInfo)

        equations(outY, fCurrentBuffer)

        while (time < endTime){
            step = normalizeStep(step, time, endTime)

            for (i in yNextBuffer.indices){
                yNextBuffer[i] = outY[i] + step*fCurrentBuffer[i]
            }

            equations(yNextBuffer, fNextBuffer)

            stepInfo.evaluationsCount++

            for (i in vectorBuffer.indices){
                vectorBuffer[i] = fNextBuffer[i] - fCurrentBuffer[i]
            }

            val errNorm = 0.5 * step * zeroSafetyNorm(vectorBuffer, outY, 1e-7)
            val q = sqrt(accuracy / errNorm)

            if(q < 1.0 && !stepInfo.isLowStepSizeReached && !stepInfo.isHighStepSizeReached){
                step = q * step / 1.1

                stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
                stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)
                stepInfo.returnsCount++
            }
            else {
                for (i in outY.indices){
                    outY[i] = yNextBuffer[i]
                }

                time += step

                stepInfo.stepsCount++

                executeStepHandlers(time, outY, stepInfo)

                step = q * step / 1.1

                stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
                stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)

                val temp = fCurrentBuffer
                fCurrentBuffer = fNextBuffer
                fNextBuffer = temp
            }
        }
        return stepInfo
    }

}