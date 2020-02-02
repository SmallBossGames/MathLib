package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.differentialEquations.IFirstOrderIntegrator
import kotlin.math.*

class StabilityControlFourthOrderIntegrator(
    override val maxEvaluations: Int,
    override val evaluations: Int,
    override val accuracy: Double
) : RungeKuttaIntegratorBase(), IFirstOrderIntegrator{

    override fun integrate(
        t0: kotlin.Double,
        y0: kotlin.DoubleArray,
        t: kotlin.Double,
        outY: kotlin.DoubleArray,
        equations: (t: kotlin.Double, inY: kotlin.DoubleArray, outY: kotlin.DoubleArray) -> kotlin.Unit
    ) {
        executeStepHandlers(t0, y0)

        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fLastBuffer = DoubleArray(y0.size)

        val yNextBuffer = DoubleArray(y0.size)

        val kMatrix = Array(5){ DoubleArray(y0.size) }

        var time = t0
        var step = t / evaluations

        equations(time, outY, fLastBuffer)

        while (time < t0 + t) {
            for (i in fLastBuffer.indices) {
                kMatrix[0][i] = step * fLastBuffer[i]
                yNextBuffer[i] = outY[i] +
                        1.0 / 3.0 * kMatrix[0][i]
            }

            equations(time, yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                kMatrix[1][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] +
                        1.0 / 6.0 * kMatrix[0][i] +
                        1.0 / 6.0 * kMatrix[1][i]
            }

            equations(time, yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                kMatrix[2][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] +
                        1.0 / 8.0 * kMatrix[0][i] +
                        3.0 / 8.0 * kMatrix[2][i]
            }

            equations(time, yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                kMatrix[3][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] +
                        1.0 / 2.0 * kMatrix[0][i] -
                        3.0 / 2.0 * kMatrix[2][i] +
                        2.0 * kMatrix[3][i]
            }

            equations(time, yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                kMatrix[4][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] +
                        1.0 / 6.0 * kMatrix[0][i] +
                        2.0 / 3.0 * kMatrix[2][i] +
                        1.0 / 6.0 * kMatrix[2][i]
            }

            val timeNext = step + time
            equations(timeNext, yNextBuffer, fCurrentBuffer)

            val q1 = findQ1(kMatrix, accuracy)
            if (q1 < 1.0) {
                step = min(t + t0 - time , q1 * step / 1.1)
                continue
            }

            val r = findR(kMatrix[0], kMatrix[1], kMatrix[2])

            step = max(step, min(q1, r)*step)

            val temp = fLastBuffer
            fLastBuffer = fCurrentBuffer
            fCurrentBuffer = temp

            time = timeNext

            for (i in yNextBuffer.indices) {
                outY[i] = yNextBuffer[i]
            }

            executeStepHandlers(time, outY)
        }
    }

    private fun findQ1(kMatrix: Array<DoubleArray>, accuracy: Double) : Double{
        var norm = 0.0
        for (i in kMatrix[0].indices){
            val d = 2*kMatrix[0][i] - 9*kMatrix[2][i] + 8*kMatrix[3][i] - kMatrix[4][i]
            norm += d*d
        }
        norm = sqrt(norm)

        //стр. 95
        return (5*accuracy.pow(5.0/4.0))/(norm/30.0)
    }

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