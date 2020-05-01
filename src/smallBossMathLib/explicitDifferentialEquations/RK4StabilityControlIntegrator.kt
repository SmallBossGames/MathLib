package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.zeroSafetyNorm
import kotlin.math.*

private const val v = 1e-7

class RK4StabilityControlIntegrator(
    val defaultStep: Double,
    val accuracy: Double
) : ExplicitIntegrator() {

    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (y: DoubleArray, outF: DoubleArray) -> Unit
    ) : IExplicitMethodStatistic {
        y0.copyInto(outY)

        val fCurrentBuffer = DoubleArray(y0.size)

        val yNextBuffer = DoubleArray(y0.size)
        val vectorBuffer1 = DoubleArray(y0.size)

        val k = Array(5){ DoubleArray(y0.size) }

        var time = t0
        var step = defaultStep

        val statistic = ExplicitMethodStatistic(
            stepsCount = 0,
            evaluationsCount = 0,
            returnsCount = 0
        )

        val state = ExplicitMethodStepState(
            isLowStepSizeReached = isLowStepSizeReached(step),
            isHighStepSizeReached = isHighStepSizeReached(step)
        )

        val endTime = t0 + t

        executeStepHandlers(t0, y0, state, statistic)

        while (time < endTime) {
            step = normalizeStep(step, time, endTime)

            equations(outY, fCurrentBuffer)
            statistic.evaluationsCount++

            for (i in fCurrentBuffer.indices) {
                k[0][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 3.0 * k[0][i]
            }

            equations(yNextBuffer, fCurrentBuffer)
            statistic.evaluationsCount++

            for (i in fCurrentBuffer.indices) {
                k[1][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 6.0 * k[0][i] + 1.0 / 6.0 * k[1][i]
            }

            equations(yNextBuffer, fCurrentBuffer)
            statistic.evaluationsCount++

            for (i in fCurrentBuffer.indices) {
                k[2][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 8.0 * k[0][i] + 3.0 / 8.0 * k[2][i]
            }

            equations(yNextBuffer, fCurrentBuffer)
            statistic.evaluationsCount++

            for (i in fCurrentBuffer.indices) {
                k[3][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 2.0 * k[0][i] - 3.0 / 2.0 * k[2][i] + 2.0 * k[3][i]
            }

            equations(yNextBuffer, fCurrentBuffer)
            statistic.evaluationsCount++

            for (i in fCurrentBuffer.indices) {
                k[4][i] = step * fCurrentBuffer[i]
            }

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = 2*k[0][i] - 9*k[2][i] + 8*k[3][i] - k[4][i]
            }

            val cNorm = accuracy.pow(5.0/4.0) / (zeroSafetyNorm(vectorBuffer1, outY, v) / 150.0)
            val q1 = cNorm.pow(1.0/4.0)

            if (q1 < 1.0 && !state.isLowStepSizeReached && !state.isHighStepSizeReached) {
                step *= q1

                state.isLowStepSizeReached = isLowStepSizeReached(step)
                state.isHighStepSizeReached = isHighStepSizeReached(step)
                statistic.returnsCount++

                continue
            }

            for (i in yNextBuffer.indices) {
                outY[i] += 1.0 / 6.0 * k[0][i] + 2.0 / 3.0 * k[3][i] + 1.0 / 6.0 * k[4][i]
            }

            time += step

            statistic.stepsCount++

            executeStepHandlers(time, outY, state, statistic)

            val q2 = cNorm.pow(1.0/5.0)
            val r = findR(k[0], k[1], k[2])

            step = max(step, min(q2, r)*step)
        }
        return statistic
    }

    /*private fun findQ1(kMatrix: Array<DoubleArray>, accuracy: Double) : Double{
        var norm = 0.0
        for (i in kMatrix[0].indices){
            val d = 2*kMatrix[0][i] - 9*kMatrix[2][i] + 8*kMatrix[3][i] - kMatrix[4][i]
            norm += d*d
        }
        norm = sqrt(norm)/150.0
        val result = accuracy.pow(5.0/4.0) / norm
        //стр. 95
        return result.pow(1.0/5.0)
    }*/

    private fun findR(k1: DoubleArray, k2: DoubleArray, k3: DoubleArray) : Double {
        require(k1.size == k2.size && k2.size == k3.size)

        var max = Double.NEGATIVE_INFINITY
        for (i in k1.indices){
            val result = abs((k3[i] - k2[i]) / (k2[i] - k1[i]))
            if(result > max)
                max = result
        }

        return 3.5/(6*max)
    }
}