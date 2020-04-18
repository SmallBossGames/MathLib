package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.Integrator
import smallBossMathLib.shared.StepInfo
import smallBossMathLib.shared.zeroSafetyNorm
import kotlin.math.*

private const val g = 1.0 / 16.0
private const val alpha2 = 1.0/3.0
private const val beta21 = 1.0/3.0
private const val beta31 = 7.0/18.0
private const val beta32 = 7.0/18.0
private const val alpha3 = 7.0/9.0
private const val p1 = 1.0/7.0
private const val p2 = 3.0/8.0
private const val p3 = 27.0/56.0

private const val v = 1e-7

class RK23StabilityControlIntegrator(val evaluations: Int, val accuracy: Double)
    : ExplicitIntegrator()
{
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (inY: DoubleArray, outY: DoubleArray) -> Unit
    ) : IExplicitMethodStepInfo {
        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fLastBuffer = DoubleArray(y0.size)
        var currentEvaluationsCount = 0

        val yNextBuffer = DoubleArray(y0.size)
        val vectorBuffer1 = DoubleArray(y0.size)

        val k1 = DoubleArray(y0.size)
        val k2 = DoubleArray(y0.size)
        val k3 = DoubleArray(y0.size)

        var time = t0
        var step = t / evaluations

        val endTime = t0 + t

        val stepInfo = ExplicitMethodStepInfo(
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step),
            stepsCount = 0,
            evaluationsCount = 0,
            returnsCount = 0
        )

        executeStepHandlers(time, outY, stepInfo)

        equations(outY, fLastBuffer)
        stepInfo.evaluationsCount++;

        while (time < endTime) {
            step = normalizeStep(step, time, endTime)

            for (i in fLastBuffer.indices) {
                k1[i] = step * fLastBuffer[i]
                yNextBuffer[i] = outY[i] + beta21 * k1[i]
            }

            equations(yNextBuffer, fCurrentBuffer)
            stepInfo.evaluationsCount++;

            for (i in fCurrentBuffer.indices) {
                k2[i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + beta31 * k1[i] + beta32 * k2[i]
            }

            equations(yNextBuffer, fCurrentBuffer)
            stepInfo.evaluationsCount++;

            for (i in fCurrentBuffer.indices) {
                k3[i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + p1 * k1[i] + p2 * k2[i] + p3 * k3[i]
            }
            
            equations(yNextBuffer, fCurrentBuffer)
            stepInfo.evaluationsCount++;

            for (i in fLastBuffer.indices){
                vectorBuffer1[i] = k2[i] - k1[i]
            }

            val q1 = ((6.0 * abs(alpha2) * accuracy / abs(1.0 - 6.0*g)) / zeroSafetyNorm(vectorBuffer1, outY, v))
                .pow(1.0/3.0)

            if (q1 < 1.0 && !stepInfo.isLowStepSizeReached && !stepInfo.isHighStepSizeReached) {
                step = q1 * step / 1.1

                stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
                stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)
                stepInfo.returnsCount++

                continue
            }

            for (i in fLastBuffer.indices){
                vectorBuffer1[i] = fCurrentBuffer[i] - fLastBuffer[i]
            }

            val q2 = ((6.0*accuracy / (abs(1.0 - 6.0*g)*step)) / zeroSafetyNorm(vectorBuffer1, outY, v))
                .pow(1.0/3.0)

            val r = findR(k1, k2, k3)

            if (q2 < 1.0 && r < 1.0 && !stepInfo.isLowStepSizeReached && !stepInfo.isHighStepSizeReached) {
                step = q2 * step / 1.1

                stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
                stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)
                stepInfo.returnsCount++

                continue
            }

            for (i in yNextBuffer.indices) {
                outY[i] = yNextBuffer[i]
            }

            time += step

            stepInfo.stepsCount++

            executeStepHandlers(time, outY, stepInfo)

            step = if(q2 < 1) min(q1, q2)*step else max(step, min(q1, min(q2, r))*step)

            stepInfo.isLowStepSizeReached = isLowStepSizeReached(step)
            stepInfo.isHighStepSizeReached = isHighStepSizeReached(step)

            val temp = fLastBuffer
            fLastBuffer = fCurrentBuffer
            fCurrentBuffer = temp

            currentEvaluationsCount++
        }
        return stepInfo
    }

    /*private fun findQ1(k1: DoubleArray, k2: DoubleArray, accuracy: Double) : Double{
        require(k1.size == k2.size)

        var norm = 0.0
        for (i in k1.indices){
            val d = (k2[i] - k1[i])
            norm += d*d
        }
        norm = sqrt(norm)

        //стр. 95
        return (((6.0 * 2.0/3.0 * accuracy) / (1.0 - 6.0 * 1.0/16.0)) / norm)
    }*/

    /*private fun findQ2(f1: DoubleArray, f2: DoubleArray, step:Double, accuracy: Double) : Double {
        require(f1.size == f2.size)

        var norm = 0.0
        for (i in f1.indices){
            val d = f1[i] - f2[i]
            norm += d*d
        }
        norm *= step

        return ((6.0*accuracy / (1.0 - 6.0 * 1.0 / 16.0)) / norm)
    }*/

    private fun findR(k1: DoubleArray, k2: DoubleArray, k3: DoubleArray) : Double {
        var max = 0.0
        for (i in k1.indices){
            max = max(abs((k3[i] - k2[i]) / (k2[i] - k1[i])), max)
        }

        return 2.0/max
    }
}